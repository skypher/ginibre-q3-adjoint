# Fundamental-minus normal form for the SU(2) Q3 target

Date: 2026-09-26

This note records one proved reduction, one conjectural strengthening with
exact bounded evidence, and five exact counterexamples to natural uniform
constructions.  It moves no theorem slot: the target `(PCP)` of
Proposition 19 in `CENTRAL_CHARACTER_Q3_SEARCH.md` remains open.

## Notation

Work in the doubled ring `R(G x G)`, `G=SU(2)`, or in the doubled fusion ring
`SU(2)_k tensor SU(2)_k`.  A doubled character is an array `F[a][b]`, the
coefficient of `V_a x V_b`.  Put

```text
S_p = V_p x 1 + 1 x V_p,        D_q = V_q x 1 - 1 x V_q,
H_q = sum_(a+b=q-1) V_a x V_b = Sym^(q-1)(C^2_x + C^2_y),   q>=1.
```

The partial character of a doubled character `R` is its column zero,

```text
G_R = sum_a R[a][0] chi_a = integral_y R(x,y) dy.
```

Proposition 19 of `CENTRAL_CHARACTER_Q3_SEARCH.md` states that the full
central-character `Q3` theorem is equivalent to `(PCP)`: `G_R` has
nonnegative coefficients for every signed word
`R = prod_i S_(p_i) prod_j D_(q_j)`.  Proposition 21 gives the same
equivalence in `SU(2)_k`.

## Lemma FM1 (fundamental-minus factorization)

For every `q>=1`, in `R(G x G)`,

```text
D_q = D_1 H_q.                                             (FM1)
```

The same identity holds in `SU(2)_k tensor SU(2)_k` for `1<=q<=k`.

**Proof.**  Put `A=V_1 x 1`, `B=1 x V_1`, and let `U_n` be the Chebyshev
recursion `U_(-1)=0`, `U_0=1`, `U_(n+1)=tU_n-U_(n-1)`.  The Clebsch--Gordan
rule `V_1V_n=V_(n+1)+V_(n-1)` (`n>=1`) gives `V_n x 1=U_n(A)` and
`1 x V_n=U_n(B)`.  In `SU(2)_k` the same rule holds for `1<=n<=k-1`, so these
formulas hold for `0<=n<=k`.  Put `h_q=sum_(a+b=q-1)U_a(A)U_b(B)`, so
`h_q=H_q` and `h_1=1`.  Splitting off `b=0` and applying the recursion to
`U_b(B)` gives

```text
h_(q+1) = U_q(A) + B h_q - h_(q-1)          (h_0=0).
```

Induct on `q`.  The case `q=1` is `(A-B)h_1=A-B`.  If `(FM1)` holds for
`q-1` and `q`, then

```text
(A-B)h_(q+1)
 = (A-B)U_q(A) + B(U_q(A)-U_q(B)) - (U_(q-1)(A)-U_(q-1)(B))
 = [A U_q(A)-U_(q-1)(A)] - [B U_q(B)-U_(q-1)(B)]
 = U_(q+1)(A)-U_(q+1)(B).
```

In `SU(2)_k` every `U_n` used has `n<=q<=k`, and every summand `V_a x V_b`
of `H_q` has `a,b<=q-1<=k-1`.  QED.

The exact replay

```text
probe_su2_fundamental_minus_cone --identity 12
probe_su2_fundamental_minus_cone --identity 12 --level 7
probe_su2_fundamental_minus_cone --identity 6  --level 3
```

checks `(FM1)` against direct Clebsch--Gordan/fusion multiplication on four
base arrays (48, 28, and 12 comparisons, no mismatch).  The previously used
identity `D_2=D_1S_1` (`SU2_V2_MINUS_FUNDAMENTAL_REDUCTION_2026_07_25.md`) is
the case `q=2`, since `H_2=S_1`.

## Corollary FM2 (fundamental-minus normal form)

Every signed word with plus labels `p_i` and `m` minus labels `q_j>=1`
satisfies

```text
prod_i S_(p_i) prod_j D_(q_j) = [prod_i S_(p_i) prod_j H_(q_j)] D_1^m.
```

Hence `(PCP)`, and therefore `Q3`, is equivalent to the statement: for every
product `F` of at most `m` factors `H_q` (`q>=2`) and any number of factors
`S_p`, the partial character of `F D_1^m` has nonnegative coefficients.  The
same holds in `SU(2)_k` with all labels at most `k`.

**Proof.**  Apply `(FM1)` to each minus factor.  Since `H_1=1`, at most `m`
factors `H_q` are nontrivial.  Conversely every such `F D_1^m` arises this
way: pad the list of `H` factors with `H_1`.  QED.

In this form all minus factors are fundamental, and `D_1^(2r)` is
pointwise nonnegative on `G x G`.  In the torus coordinates of
Proposition 13, `D_1=S_1(u)S_1(v)`, and every `H_q` has the positive feature
expansion `H_q=sum_(j in W_q) chi_(j-1)(u)chi_(j-1)(v)`.

## Conjecture FM3 (enlarged fundamental-minus cone)

Drop the bound on the number of `H` factors: for every finite product `F` of
factors from `{S_p:p>=1}` together with `{H_q:q>=2}`, and every `m>=0`,
`G_(F D_1^m)` has nonnegative coefficients (in `SU(2)_k`, labels at most `k`).

By Corollary FM2, FM3 implies `(PCP)` and hence `Q3`.  FM3 is strictly
stronger as a statement, because `H_q` is not itself of the form `S_p` or a
product of such factors.

Exact bounded evidence (`boost::multiprecision::cpp_int`):

```text
--cone 4 5 6              792 products,   5,544 checks, 0 failures
--cone 6 6 8           12,376 products, 111,384 checks, 0 failures
--cone 8 7 4 --onlyH    3,432 products,  17,160 checks, 0 failures
--cone 6 6 8 --level k  k=2,...,8,10: 0 failures at every level
--cone 3 3 4 --control  adjoining X_1=V_1 x V_1+1: 51 failures (negative control)
```

Here `--cone L N M` enumerates every product of at most `N` generators with
labels at most `L`, and every `0<=m<=M`.  This is bounded discovery
evidence, not a proof.

## Exact kills of natural uniform constructions

Each item was a candidate mechanism for proving FM3 or `(PCP)` uniformly.
Each is refuted by the displayed exact counterexample.

- **K1 (super Schur--Weyl, termwise).**  Write
  `D_1^m=sum_lambda f^lambda sch S_lambda(C^(2|2))` with
  `C^(2|2)=C^2_x + Pi C^2_y`.  Positivity of each
  `G_(F sch S_lambda)` separately is false: `F=S_1`, `lambda=(1,1,1)` gives
  coefficient `-1` at `V_0`.  (`probe_su2_fm_hookwise.py 3 3 5`: 321 failures
  among 1,008 checks.)
- **K2 (Sp(4)-irreducible plus factors).**  `H_q=Sym^(q-1)C^4` are
  restrictions of `Sp(4)`-irreducibles, but not every `Sp(4)`-irreducible is
  admissible: `Res Lambda^2_0 C^4=V_1 x V_1+1` with `m=2` gives coefficient
  `-1` at `V_2` (`probe_su2_fm_sp4_cone.py 3 3 6`).
- **K3 (GL(4)-irreducible plus factors).**  `Lambda^2 C^4=V_1 x V_1+2` fails
  at `m=3` (coefficient `-1` at `V_3`); so do `S_(2,1)`, `S_(3,1)`, and
  `S_(2,1,1)`.  Moreover `S_(2,2)C^4` is admissible through `m=7`, but
  `S_(2,2)C^4 S_1` fails at `m=3`.  Hence the full admissible cone is not
  closed under multiplication by `S_1`, and any induction on factors must use
  a proper subcone (`probe_su2_fm_gl4_cone.py 4 7`).
- **K4 (walk-monotone ballot transfer).**  In the quadrant-walk form
  `[G_(F D_1^m)]_p = sum_(a,b) (-1)^b F_(ab) W_m((a,b)->(p,0))`, both local
  lemmas of the natural ballot argument are false.  Monotonicity
  `W_m((a+1,b-1)->(p,0))>=W_m((a,b)->(p,0))` fails at `m=2`, `(1,1)->(0,0)`
  (`1<2`).  The per-antidiagonal ballot condition is not preserved by
  multiplication, because `S_1^2` violates it although `S_1^2` is admissible
  (`probe_su2_fm_ballot_lemmas.py 10 4 4`).
- **K5 (Koszul complex of the abelian block `n=Hom(C^2_x,C^2_y)` of gl(4)).**
  For `F=tensor_j Sym^(k_j)C^4` the operators of `n` commute, giving a complex
  `F tensor (C^2_x -> C^2_y)^(tensor m)` with `d^2=0` whose Euler
  characteristic is `F D_1^m`.  For a single factor, the `SU(2)_y`-invariant
  cohomology has no odd part in the tested cases
  `(k,m)=(1,1),(2,1),(1,2),(2,2),(3,2),(2,3)`.  For
  `F=Sym^1 tensor Sym^1=S_1^2`, however, `H^1` contains `V_2` at `m=2` and
  `V_1+2V_3` at `m=3` (`probe_su2_fm_nkoszul.py 2 1 1`; ranks mod
  `2^31-1`).  The failure at `m=3` is forced for a whole class.  Every
  `GL(4)`-isotypic component of `F` is a `gl(4)`-submodule, so any
  differential built from the `n`-action preserves the splitting
  `S_1^2=Sym^2C^4+Lambda^2C^4`.  The Euler characteristic of the
  `Lambda^2` summand at `m=3` is `2V_1-V_3`, which is checked directly in
  the doubled ring.  Hence its odd invariant cohomology is nonzero for every
  differential that preserves the `GL(4)`-isotypic splitting.  The full
  Euler characteristic is `3V_1+2V_3+V_5`.

A bounded pattern was also observed.  For symmetric kernels supported on one
antidiagonal `a+b=n` (`n<=8`, `m<=9`), admissibility coincides with
nonnegativity of the partial sums of `(-1)^b F_(n-b,b)` from `b=0`
(`probe_su2_fm_antidiagonal.py 8 9`).  K4 shows that this rule does not
extend to several antidiagonals.

## Active candidate and first unproved formula

The live target is FM3, or equivalently the bounded-`H` form of FM2.  The
next construction must be an odd `SU(2)_x`-equivariant differential on the
invariant complex that mixes `GL(4)`-isotypic components of the plus part and
also acts on the non-`n`-module factors `S_p`, `p>=2`.  The smallest
instance forced by K5 is `F=S_1^2`, `m=3`.  There the differential must pair
the `Lambda^2C^4` summand with the `Sym^2C^4` summand, so that the odd
invariant cohomology vanishes and the even part has character
`3V_1+2V_3+V_5`.

## Gate receipt

- Selected package: uniform construction for `(PCP)` (Proposition 19).
- Locked criterion: `(PCP)` for every signed word, all levels.
- First unresolved statement: FM3, or `(PCP)` itself.
- Work classification: FM1/FM2 reduction (proved); K1--K5 direct
  construction attempts on the locked statement, each killed by an exact
  counterexample.
- Hard-target meter: unchanged.  No theorem slot moved.
- Non-moving auxiliary batches: 0 (no prerequisite detour was opened).

## Second pass on Conjecture FM3

### Falsification attempts

The randomized exact tester `stress_su2_fm3_random` (OpenMP,
`boost::multiprecision::cpp_int`) samples random products of `S_p`, `H_q`
and checks every coefficient of `G_(F D_1^m)` for all `m<=MMAX`:

```text
SAMPLES NMAX LMAX MMAX  options                  checks       failures
 2000    6    8    8                               478,944     0
20000    8   10   12   --hprob 1.0              10,653,552     0
20000   10    6   14   --hprob 0.8 --seed 7     11,253,255     0
20000   16    6    5   --hprob 1.0 --seed 11     4,908,108     0
20000   14   12    4   --hprob 0.7 --seed 13     5,757,000     0
20000   20    3    6   --hprob 0.6 --seed 17     4,488,925     0
 2000    6    6    6   --control                   296,877     23,226
```

The fourth row is the enlarged regime: products of up to 16 factors `H_q`
against at most five fundamental minus factors.  The last row is the negative
control, which adjoins `V_1 x V_1+1` to half of the samples.  This is bounded
evidence.  Exact cancellations are very deep: the signed value is routinely
about `10^-13` of its positive part, so no ratio-based margin is informative.

### Lemma FM4 (scalar form)

FM3 is equivalent to

```text
ell_r(F) := [V_0 x V_0] F D_1^(2r) >= 0     for every F in the cone and r>=0.
```

The same holds in `SU(2)_k`.

**Proof.**  The scalar value is the `V_0` coefficient of `G_(F D_1^(2r))`.
Conversely, fix a target `p>=1`.  Every generator is symmetric under
exchanging the two tensor factors.  If `m` is even, then

```text
[V_0 x V_0] S_p F D_1^m = 2 [V_p] G_(F D_1^m),
```

and `S_p F` is in the cone.  If `m` is odd, `F D_1^m` is antisymmetric, so
`[V_0 x V_0] D_p F D_1^m = 2 [V_p] G_(F D_1^m)`.  By `(FM1)`,
`D_p F D_1^m = H_p F D_1^(m+1)`, and `H_p F` is in the cone (for `p=1`, read
`H_1=1`).  For `p=0` and odd `m` the coefficient vanishes by
antisymmetry.  The finite-level argument uses the fusion trace in the same
way (Proposition 21).  QED.

### Lemma FM5 (Kostant multiplet form of FM1)

Identify `K=SU(2) x SU(2)=Sp(2) x Sp(2)` with the Levi-type subgroup of
`Sp(4)` stabilizing `C^4=C^2_x+C^2_y`.  The two groups share a maximal torus,
and the `C^4` weights are `z^(+-1)`, `w^(+-1)`.  For every `Sp(4)`-dominant
`lambda=(lambda_1>=lambda_2>=0)`,

```text
D_1 Res V_lambda^(Sp4) = V_(lambda_1+1) x V_(lambda_2)
                         - V_(lambda_2) x V_(lambda_1+1).        (FM5)
```

**Proof.**  With `rho=(2,1)`, the `C_2` Weyl alternant is
`A_mu = det[x_j^(mu_i)-x_j^(-mu_i)]`, where `(x_1,x_2)=(z,w)`.  Direct
factorization gives

```text
A_rho = (z-z^-1)(w-w^-1)[(z+z^-1)-(w+w^-1)] = A^K_(rho_K) D_1.
```

Hence `D_1 chi_lambda = A_(lambda+rho)/A^K_(rho_K)`.  Expanding the `2 x 2`
determinant `A_(lambda+rho)`, with `lambda+rho=(lambda_1+2,lambda_2+1)`, and
dividing by `(z-z^-1)(w-w^-1)` gives the two `K`-characters.  QED.

`(FM1)` is the case `lambda=(q-1,0)`, since `Res Sym^(q-1)C^4=H_q`.
`(FM5)` is the equal-rank Dirac-index (Gross--Kostant--Ramond--Sternberg)
multiplet for `(Sp(4),K)`.  Its spinor modules are `S^+=C^2 x 1` and
`S^-=1 x C^2`, so `D_1=S^+-S^-`.  Since `R(K)` is a polynomial ring, `(FM5)`
shows that every antisymmetric doubled character is `D_1` times the
restriction of a unique virtual `Sp(4)` character.

**Corollary FM6.**  `D_1^2=|A_rho/A^K_(rho_K)|^2` on the torus, so for every
exchange-symmetric `F`

```text
ell_1(F) = 2 <F^, 1>_(Sp(4)),     ell_r(F) = 2 <F^ (D_1^2)^(r-1), 1>_(Sp(4)),
```

where `F^` is the unique virtual `Sp(4)` character with `Res F^=F`.  Also

```text
S_1^ = V_(1,0),   S_p^ = V_(p,0) - V_(p-1,1) + V_(p-2,0)  (p>=2),   equivalently
S_p  = 2(H_(p+1)+H_(p-1)) - S_1 H_p   (H_0 := 0),
D_1^2^ = V_(2,0) - 3V_(1,1) + 5.
```

The formula for `S_p^` follows from `(FM5)`: multiply both sides by `D_1`,
compare with `D_1S_p`, and cancel `D_1` in the domain `R(K)`.  The
`sp4_decomposition` probe confirms it for `p<=7`.

**Corollary FM7 (two-minus sector).**  For a word with plus labels `p_i` and
two minus labels `q,q'>=1`,

```text
J = 2 < prod_i S_(p_i)^ tensor Sym^(q-1)C^4 , Sym^(q'-1)C^4 >_(Sp(4)).
```

Thus the complete two-minus sector of `Q3` (any number of plus factors) is
the statement that the virtual `Sp(4)` character `prod_i S_(p_i)^` pairs
nonnegatively with every `Sym^a tensor Sym^b`.  The repository proves this
sector in general only through seven factors (Corollary 23A9ZZ10), plus the
special arbitrary-length suffix cones of Propositions 24A--24B2.

### Lemma FM8 (Andreief closed form; integrality is necessary)

Put `alpha=(theta+phi)/2`, `beta=(theta-phi)/2` for the two `SU(2)` angles,
and `tau=cos^2 alpha`, `tau'=cos^2 beta`.  Then:

- `sin theta sin phi = tau' - tau`;
- `D_1 = -4 sin alpha sin beta`;
- `F` splits into parts even-even and odd-odd in `(cos alpha, cos beta)`;
- the even-even part `F_ee(tau,tau')` is a symmetric polynomial.

Write `F_ee = sum_lambda c_lambda(F) s_lambda(tau,tau')` in two-variable Schur
polynomials.  With Andreief's identity and the Beta moments
`int tau^k (1-tau)^r d(arcsine) ∝ (1/2)_k/(r+1)_k`,

```text
ell_r(F)/ell_r(1) =
  [ sum_lambda c_lambda (lambda_1-lambda_2+1) (1/2)_(lambda_1+1) (1/2)_(lambda_2)
        / ((r+1)_(lambda_1+2) (r+1)_(lambda_2+1)) ]
  / [ (1/2) / ((r+1)_2 (r+1)_1) ].
```

The probe `probe_su2_fm3_closed_form.py 5` checks this identity exactly
against direct doubled-ring multiplication for `r<=5` and six products.
Examples:

```text
H_3:        2r(r-1)/((r+2)(r+3))
S_2:        6r(r+1)/((r+2)(r+3))
H_5:        3r(r-1)^2(r-2)/((r+2)(r+3)^2(r+4))
S_1^2 H_3:  24(r^2-3r+12)/((r+2)(r+3)^2(r+4))
```

The natural real interpolation in `r` is therefore negative on `(0,1)` for
`H_3` and on `(1,2)` for `H_5`.  Any proof of FM3 must use the integrality of
`r`.  This excludes, for instance:

- Schur positivity in `tau`;
- positivity of Selberg or Kadell averages at real parameters;
- any Cauchy--Schwarz argument over a real family of weights.

Every `s_lambda(tau,tau')` has strictly positive average for all real
`r>-1/2`, so `F_ee` is never Schur-positive in `tau` when a zero at
`r=1` occurs.

Two equivalent ensemble pictures come out of the same change of variables.
In the `SU(2)` cosines `c=cos theta`, `ell_r` is the average over the
two-point Jacobi `beta`-ensemble with `beta=2r` and semicircle weight.  In
`tau`, it is the `beta=2` ensemble with weight
`tau^(-1/2)(1-tau)^(r-1/2)`.  At `r=1` both reduce to `Sp(4)` Haar measure.

### Further bounded observations and kills

- **Sp(2n) analogue.**  The same statement with `n` copies holds in every
  tested case: `n=3` (labels `<=3`, three generators, `r<=2`; labels `<=4`,
  four generators, `r=1`) and `n=4` (labels `<=2`, three generators, `r=1`).
  Here the cone is generated by `sum_i V_p(x_i)` and
  `Res Sym^(q-1)C^(2n)`, and the weight is
  `prod_(i<j)(chi_1(x_i)-chi_1(x_j))^(2r)`.  The script is
  `probe_su2_fm3_ncopy.py`.
- **K6.**  Splitting `D_1^(2r)=sum_lambda f^lambda s_lambda(x,-y)` by
  `GL(4)` Schur--Weyl and requiring each piece to be positive fails:
  `lambda=(1,1)`, `F=H_3` gives `-1` (`probe_su2_fm3_gl4_twisted.py 6 4 4`).
- **K7.**  The `Sp(4)` multiplicities of products `prod S_p^` obey no sign
  rule by the parity of `lambda_2`.  Among 34 products with labels `<=4`
  and at most three factors, 15 violate each candidate parity rule.  Only
  the trivial multiplicity stays nonnegative
  (`probe_su2_fm3_sp4_products.py 4 3`).

### Current first unproved statement

The smallest open instance is the `r=1`, `H`-free part of FM3:

```text
< prod_i S_(p_i)^ , 1 >_(Sp(4)) >= 0,
S_p^ = V_(p,0) - V_(p-1,1) + V_(p-2,0).
```

Equivalently, `Q3` holds for words with minus multiset `{1,1}` and
arbitrary plus labels.  It is not covered by the repository's
`V_1/V_2` sector theorem, which requires every label to be at most two.

## Third pass: the two-fundamental-minus sector

### Lemma FM9 (exponential Heine form)

For a multiset `alpha` of plus labels (`alpha_p` copies of `p`) and `2r`
fundamental minus factors, put

```text
F_s = exp(sum_p s_p chi_p),    f_a(s) = [chi_a] F_s,    M_k(s) = <F_s chi_1^k>.
```

Then

```text
J(alpha; 1^(2r)) = [s^alpha/alpha!] sum_(k=0)^(2r) (-1)^k binom(2r,k) M_k M_(2r-k).
```

For `r=1` this is

```text
J(alpha; 1,1) = 2 [s^alpha/alpha!] ( f_0 (f_0+f_2) - f_1^2 ).
```

**Proof.**  Since `exp(sum_p s_p(chi_p(x)+chi_p(y))) = F_s(x)F_s(y)`, the
coefficient of `s^alpha/alpha!` in
`int int F_s(x)F_s(y)(chi_1(x)-chi_1(y))^(2r)` is
`int int prod_p (chi_p(x)+chi_p(y))^(alpha_p) D_1^(2r)`.  Expanding
`(chi_1(x)-chi_1(y))^(2r)` gives the first formula.  For `r=1`, use
`M_0=f_0`, `M_1=f_1`, `M_2=f_0+f_2` (because `chi_1^2=chi_2+chi_0`) and the
symmetry of the sum.  QED.

Equivalently, by Andreief's identity with the functions `(1,chi_1)`,
`J/2 = [s^alpha/alpha!] det[<F_s chi_j chi_k>]_(j,k<2)`, the Sp(4) Heine
determinant (corrected by FM-CHK; at `s=0`, `J=2` and the determinant is `1`).
The `Sp(2n)` version of the same sector is the `n x n` Hankel determinant
`det[M_(j+k)]_(j,k<n)=det[<F_s chi_j chi_k>]_(j,k<n)`.  The Sp(2n) analogue
of the `r=1` statement is therefore coefficientwise positivity of all
Hankel determinants of the tilted semicircle moments `M_k(s)`.

### Lemma FM10 (twisted normal form) and Corollary FM11 (a proved family)

Translating `y` by the central element `-1` preserves Haar measure and maps
`chi_a(y)` to `(-1)^a chi_a(y)`.  For every signed word it therefore
exchanges the sign type of every odd label.  After this twist, apply
`(FM1)` to every minus factor and `(FM5)`/Corollary FM6.  Every word then
has the form

```text
J = < Res X , D_1^(2r') >_K,    X in R(Sp(4)),
2r' = #(odd plus labels) + #(even minus labels),
```

where `X` is a product of the following factors:

- `S_p^` for every even plus label `p>=2` (virtual);
- `S_q^` for every odd minus label `q>=3` (virtual);
- `C^4` for every minus label `1`;
- `Sym^(p-1)C^4` for every odd plus label `p`;
- `Sym^(q-1)C^4` for every even minus label `q`.

In the untwisted picture instead, `2r=#minus`, the virtual factors are the
plus labels `>=2`, and every minus label contributes `Sym^(q-1)C^4`.

**Corollary FM11.**  `Q3` holds for every signed word in either of the
following classes:

1. There is no even plus label `>=2`, every odd minus label equals `1`, and
   `#(odd plus labels) + #(even minus labels) = 2`.
2. There are exactly two minus factors and every plus label equals `1`.

**Proof.**  In the respective picture `r'=1` (resp. `r=1`) and `X` is a
genuine `Sp(4)`-module.  By Corollary FM6, `J = 2 dim X^(Sp(4)) >= 0`.  QED.

Class 1 contains words of unbounded length, for example `[3^+,5^+,1^-^(2k)]`
or `[p^+ (odd), q^- (even), 1^-^(odd count)]`.  These words are not covered
by the `V_1/V_2` sector (labels `<=2`), the parity-separated chamber
(Proposition 20), or the seven-factor two-minus theorem.  The general
two-fundamental-minus sector remains open.  Its obstruction is exactly the
virtual factors `S_p^`, `p>=2`.

The same twist also shows that FM0 with all plus labels odd is equivalent to
an `H`-only instance of FM3.  Here FM0 denotes the part of FM3 with
fundamental minus factors and plus factors `S_p` only.

### Kills recorded in this pass

- **Monotonicity in n.**  The Sp(2n) values `E_n` are not monotone in `n`:
  `(1,3,4,4)` gives `E_1,E_2,E_3 = 2,7,6`.  There are three such words among
  125 (`probe_su2_fm3_sp2n_mono.py 4 5 3`).
- **Single-label total positivity.**  `exp(sA_1)` is coefficientwise `TP_2`
  (Karlin--McGregor).  For `p>=2`, `exp(sA_p)` is not: for `p=2,3,4`,
  respectively 72, 124 and 84 of 441 minors have a negative coefficient
  (`probe_su2_fm3_tp2.py`).  The restriction to the leading two columns also
  fails, for parity reasons.  The second compound of `exp(L)` is exactly the
  matrix of multiplication by `exp(sum s_p S_p^)` in the `Sp(4)`-irreducible
  basis, so this route is equivalent to the original question.
- **Row domination.**  For every tested product, the `r=1` statement is the
  family `y_(p-1,1) <= y_(p,0)+y_(p-2,0)` on `Sp(4)` multiplicities, and it
  holds on all 125 products.  Its extension
  `y_(a,b) <= y_(a+1,b-1)+y_(a-1,b-1)` to higher rows fails for both parities
  of `b` (`probe_su2_fm3_row_domination.py 5 4`).
- **Scalar positivity of single GL(4) irreducibles.**  This fails:
  `ell_r(Lambda^2 C^4)` for `r=0,...,6` is `2,2,4,0,-168,-2772,-37752`, and
  `S_(2,1,1)`, `S_(4,1,1)`, `S_(3,2,1)`, `S_(2,2,1,1)` also go negative
  (`probe_su2_fm3_gl4_scalar.py 6 6`).  The reason is that
  `D_1^(2r)/int D_1^(2r)` converges weakly to the average of the point
  masses at `epsilon=(1,-1)` and `(-1,1)`.  Hence
  `ell_r(F)/ell_r(1) -> F(epsilon)`, and
  `s_(1,1)(1,1,-1,-1) = -2`.  Every generator of the cone satisfies
  `F(epsilon) >= 0` (`S_p(epsilon)=2(p+1)[p even]` and
  `H_q(epsilon)=((q+1)/2)[q odd]`), which is necessary for FM3.
- **Jacobi continued fraction.**  The coefficients `beta_n,gamma_n` of the
  normalized tilted moments are not coefficientwise nonnegative, already for
  `p=1` (`probe_su2_fm3_jfraction.py`).  Sokal-type continued-fraction
  positivity does not apply in this normalization.

## Fourth pass: the H-only sector and a finite Sp(4) form of FM3

### Lemma FM12 (walk-weighted Sp(4) formula)

For every virtual `Sp(4)`-module `M` with multiplicities `m_mu(M)` and every
`r>=1` (the formula uses `N_(2r-1)`),

```text
ell_r(Res M) = 2 sum_mu (-1)^(mu_2) N_(2r-1)(mu_1+1, mu_2) m_mu(M),
```

where `N_n(a,b)` is the number of length-`n` walks in `N^2` with unit
steps `+-e_1,+-e_2` from `(0,0)` to `(a,b)`.  Explicitly,
`N_n(a,b) = sum_k binom(n,k) c(a,n-k) c(b,k)`, with ballot numbers `c`.

**Proof.**  Since `D_1` is real, `ell_r(Res M) = <D_1 Res M, D_1^(2r-1)>_K`.
By `(FM5)`, `D_1 Res V_mu = V_(mu_1+1) x V_(mu_2) - V_(mu_2) x V_(mu_1+1)`.
Expanding `(V_1 x 1 - 1 x V_1)^n` shows
`[V_a x V_b] D_1^n = (-1)^b N_n(a,b)`.  For odd `n`, `a+b` is odd, so the
two terms contribute equally.  QED.

Since `N_(2r-1)(a,b)=0` for `a+b>2r-1`, `ell_r` depends only on the finitely
many multiplicities `m_mu`, `|mu|<=2r-2`.  The weight vectors are:

```text
r=2: 5 m_(0,0) - 3 m_(1,1) + m_(2,0)
r=3: 35 m_(0,0) - 35 m_(1,1) + 14 m_(2,0) + 10 m_(2,2) - 5 m_(3,1) + m_(4,0)
r=4: 294 m_00 - 378 m_11 + 168 m_20 + 189 m_22 - 105 m_31 - 35 m_33
     + 27 m_40 + 21 m_42 - 7 m_51 + m_60
```

FM3 at level `r` is therefore the statement that this finite signed
functional is nonnegative on every product of the factors `S_p^` and
`Sym^q C^4`.

### Corollary FM13 (H-only, r=2, invariant-count form)

For `M = tensor_j Sym^(k_j)C^4`, let `I(k)` be the dimension of the
`Sp(4)`-invariants of `M`.  Then

```text
ell_2(M)/2 = 4 I(k + (2)) - 3 I(k + (1,1)) + 8 I(k),
```

where `k+(2)` appends one factor `Sym^2` and `k+(1,1)` appends two factors
`C^4`.

**Proof.**  `C^4 tensor C^4 = Sym^2 + Lambda^2_0 + 1`.  Hence
`I(k+(1,1)) = m_00+m_11+m_20` and `I(k+(2)) = m_20`.  Substitute into
Lemma FM12.  QED.

A direct check on 110 multisets agrees exactly.  The inequality depends on
the rank of `Sp(4)`.  For `k=1^6` the `Sp(4)` counts
`(I(k), I(k+(1,1)), I(k+(2))) = (14, 84, 40)` give `+20`.  The free
multigraph counts `(15, 105, 45)`, which are the stable large-rank values
without the `6 x 6` Pfaffian relations, give `-15`.  Any proof of FM3 must
therefore use the rank-two structure, for example 3-crossing-free standard
monomials or `C_2` Weyl-chamber walks.  Arguments valid uniformly in the
rank cannot work.

### Lemma FM14 (vanishing half of the support rule)

Sort `k_1 >= k_2 >= ...` and put `s = sum_(j>=2) k_j`.  If `r>=1` and
`k_1 >= s+2r`, then `ell_r(tensor_j Sym^(k_j)C^4) = 0`.  (The restriction
`r>=1` is necessary: `k=(1,1)` has `k_1=s` but `ell_0(C^4 tensor C^4)=1`.)

**Proof.**  If `V_mu` is a constituent of `Sym^(k_1) tensor N`, with `N` the
product of the other factors, then `Sym^(k_1)` is a constituent of
`V_mu tensor N`, since every module involved is self-dual.  So `k_1` is the
first coordinate of a weight of `V_mu tensor N`, which is at most
`mu_1+s`.  Hence `mu_1 >= 2r`, and `N_(2r-1)(mu_1+1,mu_2)=0`.  QED.

The exact probe `probe_su2_fm3_honly_support.py 12 8 4` checks 151
multisets for `r<=4`.  It finds that `ell_r>0` in every case with
`k_1 < s+2r` and the right parity, and no negative value.  The positive half
of the rule is the H-only case of FM3.

### Rational generating functions for few factors

With `G_J(t) = prod_(i<j)(1-t_i t_j)^(-1)`, the constant-term computation
(`probe_su2_fm3_qr_rational.py`) gives, for `J<=3` factors:

```text
Q_1 = 2 G_J,   Q_2 = 2 G_J (5 + h_2 - 3 e_2),
Q_3 = 2 G_J (35 + 14 p_2 - 21 e_2 + (t_1-t_2)^4)   (J=2).
```

In general `Q_r = 2 G_J sum_mu (-1)^(mu_2) N_(2r-1)(mu_1+1,mu_2) s_mu(t)`,
from Littlewood's symplectic Cauchy identity.  Three-row universal
characters vanish for `Sp(4)`, so the formula is exact for `J<=3`.  For
`J=3`, `G_3` is the indicator of the parity-restricted triangle cone `T`.
Monomial positivity at `r=2` reduces to

```text
5 1_T(k) + sum_i 1_T(k-2e_i) >= 2 sum_(i<j) 1_T(k-e_i-e_j),
```

which holds by a short case check (interior, boundary facet, exterior).
These few-factor cases lie inside the already proved `<=6`-factor range.
They are recorded for the mechanism, not as new theorem coverage.

## Fifth pass: invariant-theory model of the r=2 H-only inequality

### Standard-monomial model

The `Sp(4)`-invariants of `C^4 tensor C^J` form `C[omega_ij]/(Pf_6)`
(first and second fundamental theorems).  Order the vertices `1<...<J`.
The multidegree-`k` piece should then have a basis of loopless multigraphs
on `[J]` with degree sequence `k` and no 3-nesting
`(a_1,b_1) ⊃ (a_2,b_2) ⊃ (a_3,b_3)`, `a_1<a_2<a_3<b_3<b_2<b_1`.  This is
the standard-monomial description coming from the Gröbner basis of the
`6 x 6` Pfaffian ideal.  The precise term-order citation (Herzog--Trung,
Jonsson--Welker) is still to be checked in context.  The exact probe
`probe_su2_fm3_nesting_model.py 10 4` confirms `I(k) = #(3-nesting-free
multigraphs)` on all 1,317 ordered degree sequences with at most five
vertices, sum at most 10 and parts at most 4.

Append new vertices at the right end.  For a 3-nesting-free `Gamma` put:

- `alpha(Gamma) <= p` if no two edges starting right of `p` are nested;
- `beta(Gamma) <= q` if no edge starts right of `q`.

The only possible new 3-nestings involve the appended edges, and a direct
case analysis gives, for `M=tensor_j Sym^(k_j)C^4`,

```text
m_(0,0) = I(k),
m_(2,0) = sum_(p<=q) #{Gamma in N(k-e_p-e_q) : alpha(Gamma)<=p},
m_(1,1) = sum_(p<q)  #{Gamma in N(k-e_p-e_q) : alpha(Gamma)<=p, beta(Gamma)<=q}.
```

The probe `probe_su2_fm3_check_x.py 10 4` checks the last two formulas
exactly on 345 degree sequences.  They imply

```text
m_(1,1) <= m_(2,0)      for every product of symmetric powers
```

(conditional on the model).  The needed inequality is
`3 m_(1,1) <= m_(2,0) + 5 m_(0,0)`.

### Data and the size of the gap

The probe `probe_su2_fm3_ratios.py 16 8 8` covers 380 multisets.  On all
of them:

- `2 m_(1,1) <= m_(2,0) + 2 m_(0,0)` holds;
- `3 m_(1,1) <= m_(2,0) + 5 m_(0,0)` holds;
- `3 m_(1,1) <= m_(2,0) + 4 m_(0,0)` fails 194 times;
- `m_(1,1) <= 2 m_(0,0)` fails 307 times.

The ratio `(m_(2,0) + 5 m_(0,0))/(3 m_(1,1))` decreases toward `1`.  The
smallest observed value is `85/83`, at `(3,3,3,3,3,1)`.  The reason is that
for large products `m_mu/m_(0,0) -> dim V_mu`, which is `1, 5, 10`.  The
target weights then give `10 + 5 - 3*5 = 0` at leading order, so positivity
is a subleading effect.  Its sign is that of the pointwise-nonnegative
weight `D_1^2` near the identity.  A combinatorial proof must therefore be
asymptotically bijective.

### Kills and trivial subfamilies

- **K8.**  The row-parity sign cone fails.  For `v_lam(M) =
  ell_r(M tensor V_lam)/2`, the sign `(-1)^(lam_2)` is violated already by
  `M = C^4 tensor C^4`, `lam=(1,1)`, where the value is `+2`
  (`probe_su2_fm3_signcone.py 4 4 2 4`: 321 violations among 2,100 values).
  No cone defined by that sign pattern is invariant under tensoring with
  symmetric powers.
- **Trivial subfamilies (proved).**  Translating by `y -> -y` sends
  `prod_(q in Q) H_q D_1^(2r)` with `|Q|=2r` and all `q` odd to the all-plus
  word `prod_q S_q`.  So those H-only values are nonnegative (Proposition 20
  in twisted form).  FM0 values whose plus labels all have even
  multiplicity are nonnegative pointwise:
  `E_(Sp4)[prod X_p^(2a_p) D_1^(2r-2)] >= 0`.  The genuinely hard H-only
  words contain even `q`, or more `H` factors than the minus count allows.

## Sixth pass: two structural identities for the r=2 weight

### Lemma FM15 (exterior-algebra form of the r=2 weight)

Let `R^5 = Lambda^2_0 C^4` be the vector representation of
`Spin(5) = Sp(4)`.  In `R(Sp(4))`,

```text
D_1^2^ = V_(2,0) - 3 V_(1,1) + 5 = sum_(j=0)^5 (-1)^(j-1) j Lambda^j(R^5)
       = + d/dt lambda_t(R^5) |_(t=-1)      (sign corrected by FM-CHK).
```

**Proof.**  `Lambda^0 = Lambda^5 = 1`, `Lambda^1 = Lambda^4 = R^5 = V_(1,1)`,
and `Lambda^2 = Lambda^3 = so(5) = sp(4) = V_(2,0)`.  Hence the sum is
`R^5 - 2V_(2,0) + 3V_(2,0) - 4R^5 + 5 = V_(2,0) - 3V_(1,1) + 5`.
Alternatively, `lambda_t(R^5) = (1+t) lambda_t(p)` with `p = V_1 x V_1`, so
the derivative at `t=-1` is `lambda_(-1)(p) = D_1^2`.  QED.

So the `r=2` statement for a module `M` reads

```text
sum_j (-1)^(j-1) j dim Hom_(Spin(5))(Lambda^j R^5, M) >= 0,
```

a degree-weighted Euler characteristic of the `Lambda^j R^5`-covariants.
For an exact Koszul-type complex `C_j`, the weighted sum
`sum_j (-1)^j j dim C_j` equals the alternating sum of the boundary
dimensions.  That is the natural place to look for a complex realizing it.

### Lemma FM16 (the r=2 integral over an SU(4) subtorus)

Embed `K = SU(2) x SU(2) ⊂ Sp(4) ⊂ SU(4)`.  The `K`-torus is the
codimension-one subtorus `(z, z^(-1), w, w^(-1))` of the `SU(4)`-torus.
Restricting the `SU(4)` roots `e_i - e_j` gives `2e_1`, `2e_2` once each and
each of `+-e_1+-e_2` twice.  Hence

```text
|Delta_SU(4)|^2 restricted to T_K = |Delta_K|^2 D_1^4,
ell_2(F) = (1/4) int_(T_K) F |Delta_SU(4)|^2    for every K-character F.
```

In the H-only case every factor `Sym^k C^4` is an `SU(4)`-module.  Writing
`K = [S(U(2) x U(2)), S(U(2) x U(2))]` and summing over the central `U(1)`
characters, Atiyah--Bott localization on `Gr_2(C^4)` identifies the
integral with

```text
ell_2(M) = sum_(n in Z) < M , sum_p (-1)^p chi_SU(4)(Gr_2(C^4), Omega^p(n)) >_SU(4).
```

Here `sum_p (-1)^p chi(Omega^p(n))` has character equal to the `W`-orbit
sum of `n omega_2`.  This localization step is sketched, not yet checked
term by term.  It places the `H`-only `r=2` inequality among equivariant
`chi_y`-genera of line bundles on the Grassmannian.  Bott- or Snow-type
vanishing for `H^q(Gr_2(C^4), Omega^p(n) tensor M^*)` is a concrete
candidate mechanism.  At `r=1` the corresponding space is
`Sp(4)/K = S^4`.  For `r>=3` no compact group with the required
root-multiplicity pattern was found.

## Seventh pass: level recursion and further stress data

### Lemma FM17 (level recursion)

Put `W_r(M) = ell_r(Res M)/2 = < M , (D_1^2^)^(r-1) >_(Sp(4))`.  Since

```text
D_1^2^ = 4 Sym^2 C^4 - 3 (C^4 tensor C^4) + 8        in R(Sp(4)),
```

every (virtual) `Sp(4)`-module `M` satisfies

```text
W_(r+1)(M) = 4 W_r(M tensor Sym^2) - 3 W_r(M tensor C^4 tensor C^4) + 8 W_r(M),
W_1(M) = dim M^(Sp(4)).
```

**Proof.**  `C^4 tensor C^4 = Sym^2 + Lambda^2_0 + 1`, so
`4Sym^2 - 3(C^4)^(tensor 2) + 8 = Sym^2 - 3Lambda^2_0 + 5 = D_1^2^`
(Lemma FM13).  Then `W_(r+1)(M) = W_r(M tensor D_1^2^)` by definition.  QED.

Thus the `H`-only part of FM3 says `T^(r-1) I >= 0` for all `r`.  Here
`I(k)` is the `Sp(4)`-invariant count of `tensor_j Sym^(k_j)C^4`, and `T`
appends one of three gadgets:

- one `Sym^2` factor (weight 4);
- two `C^4` factors (weight -3);
- nothing (weight 8).

In the 3-nesting-free multigraph model, `T^(r-1)I(k)` is therefore a signed
count of multigraphs on `k` plus appended gadget vertices.  `T1` is the case
`r=2`.

### Stress data

`stress_su2_fm3_random` in the `H`-heavy regime:

- `100000 14 8 16 --hprob 1.0 --seed 101`: 94,222,687 checks, 0 failures;
- `100000 12 6 20 --hprob 0.85 --seed 202`: 99,727,047 checks, 0 failures.

Every coefficient of `G_(F D_1^m)` was checked for `m<=16` and `m<=20`
respectively, which covers every scalar level `r<=10` through targets.

## Agent returns (gpt-6-luna, reviewed by the main agent)

- **K9 (FM-T1-KZ, luna_max_saturn; accepted after re-check).**  Consider any
  complex whose terms are `Hom_(Sp(4))(Lambda^(d+1)R^5, M)^((d+1))`, so that
  its Euler characteristic is `5m_00-3m_11+m_20` by Lemma FM15.  If its
  differentials are natural in the `Sp(4)`-module `M`, then Yoneda leaves
  only `Hom(Lambda^2R^5, Lambda^3R^5) = C` (Hodge star) as a possible nonzero
  differential.  Hence `dim H_3 = 4 m_11`, which is already `4` at
  `k=(1,1)`.  For `(C^4)^(tensor 4)` the chain dimensions are
  `(5,12,18,20,15)`.  A positivity proof by odd vanishing must therefore use
  maps that depend on the symmetric-power structure of `M`.
- **FM-T2 (luna_max_venus).**  The two-plus-label pairings
  `<S_p^ S_q^,1> = 1 (p=q=1), 3 (p=q>=2), 1 (|p-q|=2), 0` are correct.
  They are four-factor words, inside the proved `<=7`-factor range, so there
  is no criterion movement.  Spot checks `(3,3,1,1,1,1) -> 11`,
  `(1,1,3,3) -> 4` and `(2,3,3) -> 6` agree.
- **New sub-target (FM-T2b).**  For words with labels in `{1,3}` in which
  all `3`s carry the same sign, `J = I(b,c,r)` (precise counts and parity
  clause as repaired by FM-CHK2 and used in Corollary FM18: all `3`s minus
  gives `b=#D_3`, `c=#S_1`, `2r=#D_3+#D_1`; all `3`s plus gives, after the
  twist, `b=#S_3`, `c=#D_1`, `2r=#S_3+#S_1`; `J=0` when that count is
  odd), where

  ```text
  I(b,c,r) = int int (x+y)^c (x-y)^(2r) (x^2+xy+y^2-2)^b dsigma(x) dsigma(y)
  ```

  and `sigma` is the semicircle law on `[-2,2]`.  The exact grid
  `b<=8`, even `c<=8`, `r<=6` has no negative value; the only zeros are
  `(1,0,0)` and `(1,0,1)`.

- **FM-CHK (luna_max_mercury), independent re-derivation of FM1--FM16.**
  - Accepted: FM1, FM2, FM4, FM5, FM7, FM8, FM10, FM11, FM13, and the first
    identity of FM16.
  - Repaired above: FM6 (`S_1^` boundary), FM9 (the Heine determinant equals
    `J/2`), FM12 (`r>=1`), FM14 (`r>=1`; false at `r=0`), FM15 (sign of the
    derivative).
  - FM10 remark: if `#(odd plus)+#(even minus)` is odd, then `J=0` by
    exchange antisymmetry.
  - The main agent re-checked each defect before editing.
- **FM-LIT (luna_max_neptune), source audit; accepted.**
  - Standard-monomial model: Jonsson--Welker, "A spherical initial ideal for
    Pfaffians", Prop. 3.1, restating Herzog--Trung, Adv. Math. 96 (1992),
    gives a term order whose initial ideal of the ideal of `6 x 6` Pfaffians
    is generated by the 3-nestings for the vertex order `1<...<n`.
    Standard monomials in multidegree `k` are therefore exactly the
    3-nesting-free loopless multigraphs, which justifies the fifth-pass model
    and hence `m_11 <= m_20`.  The Herzog--Trung variable ranking itself was
    not re-verified from the 1992 original.  Jonsson--Welker's own order
    (Thm. 2.1) gives the 3-crossing-free basis instead.
  - No audited source (Bakry--Echerbault, Ginibre, Bakry--Huet, Dunlop 1976,
    Herbst 2022, Trimeche, Opdam) settles any T1/T2 case.
  - One slip in the packet: `<S_1^2,1>_(Sp(4)) = 1`, not `2`.
- **T2 in walk form (assigned as FM-T2W).**  With `Phi_ab` the number of
  lattice walks from `(0,0)` in which step `i` moves exactly one coordinate by
  the `V_(p_i)` Clebsch--Gordan rule, T2 is exactly
  `Phi_11 <= Phi_00 + Phi_20`.
- **FM-T1-GR (luna_max_jupiter); accepted.**
  - Lemma FM16's Grassmannian identity holds with explicit conventions:
    `O(1) = det S^*` on `X = Gr_2(C^4)`,
    `E_n = sum_p (-1)^p chi_SU(4)(X, Omega^p(n))`, character of
    `E_n = sum_(|I|=2) t_I^(-n)`, and
    `ell_2(Res M) = sum_(n in Z) <M, E_n>_(SU(4))`, a finite sum.
  - Proof: Atiyah--Bott holomorphic Lefschetz (M. F. Atiyah and R. Bott,
    "A Lefschetz fixed point formula for elliptic complexes: II.
    Applications", Ann. of Math. 88 (1968), 451--491, Thm. 4.12) at a generic torus
    element, then the `u`-constant-term argument that restricts to `T_K`.
  - Numerical checks: `ell_2 = 10, 2, 6, 6` for `M = 1`, `Sym^2`,
    `(C^4)^(tensor 2)` and `(Sym^2)^(tensor 2)`.
  - T1 for at most three factors (no new coverage).
- **K10 (FM-T1-IC, luna_max_saturn).**  No simple domination family is
  invariant under the level operator `T`; `F_T` fails via an explicit
  degree-two function.  Grid data: `W_r >= 0` on 272 partitions with
  `|k| <= 12`, `r <= 5`.
- **K11 (luna_max_jupiter).**  Termwise Bott vanishing fails inside the
  cone.  For `M = V^(tensor 4)` and `(n,p,q) = (2,1,0)`,
  `H^0(Omega^1(2)) = S_(2,1,1)V` has multiplicity 3, contributing `-3`,
  which cancels within `<M, E_2> = 0`.
- **K12 (main agent).**  The `u`-refinement is not positive: `<M, E_n> >= 0`
  for each `n` fails at `k = (3,1^7)`, with coefficient `-1` at `n = +-3`
  (`probe_su2_fm3_nwise.py`).
- **New assignments.**
  - FM-T1-SX: four minus factors (arbitrary labels) plus a fundamental plus
    suffix of any length.  Exact screen: 822 cases with `q_j <= 6`,
    `a <= 12`, none negative.
  - FM-T2-123: T2 on plus labels `{1,2,3}`.  Grid `a <= 6`, `b <= 5`,
    `c <= 5`: none negative.
- **FM-T2W (luna_max_neptune).**  T2 for any three plus labels is proved
  (margin `2 tau_0 + tau_2 + E - 2N >= 0`).  These are five-factor words,
  inside the proved `<=7`-factor range, so there is no criterion movement.
  T2 is open only from six plus labels on.
- **K13.**  At four labels, grouping the walk count by subset size `|A|`
  cannot be termwise positive.  For `(1,2,2,3)` the per-size margins are
  `(2, 0, -3, -2, 7)`, with total `4`.
- **Scope note.**  Prop. 24A of the central search already proves two
  arbitrary minus labels with a fundamental plus suffix, via Karlin--McGregor
  in the wedge basis `e_u ^ e_v`.  Under FM5 that basis is the `Sp(4)`
  irreducible basis.  New assignment FM-T2-SX asks for two arbitrary PLUS
  labels with a fundamental plus suffix, `<S_p^ S_q^ (C^4)^(tensor a), 1> >= 0`.
  Exact screen: `2 <= p <= q <= 8`, `a <= 14`, 212 cases, no negative value.

### Corollary FM18 (words in {1,3} with a single 3; FM-T2b, luna_max_venus; accepted)

**Moment formula.**  Let `x,y` be independent semicircle variables on
`[-2,2]` and put `u=x+y`, `v=x-y`.  Then

```text
M(m,n) := E[u^(2m) v^(2n)]
        = 2 (2m)!(2m+1)!(2n)!(2n+1)! / ((m!)^2 (n!)^2 (m+n+1)! (m+n+2)!).
```

Venus derived this from semicircle integration by parts
`int ((4-x^2) f' - 3x f) dsigma = 0` together with the `u<->v` symmetry.
An independent derivation by the main agent: in half-angle coordinates
`u = 4 cos(alpha) cos(beta)` and `v = -4 sin(alpha) sin(beta)`, so the
integrand factorizes against the weight `(sin^2 alpha - sin^2 beta)^2`, and
Andreief reduces `M` to Beta integrals.  The formula agrees exactly with
direct expansion for `m,n <= 8`.

**b=1.**  Since `H_3 = (3u^2+v^2-8)/4`,

```text
I(1,2m,r) = M(m,r) (10m^2 - 4mr + 2r^2 + 14m - 2r) / ((m+r+2)(m+r+3)),
10m^2 - 4mr + 2r^2 + 14m - 2r = 2(r-m-1/2)^2 + 8m^2 + 12m - 1/2.
```

The right side is positive for `m>=1` and equals `2r(r-1) >= 0` for
`m=0`.  This agrees exactly with direct computation for `m,r <= 7`.

**Corollary FM18.**  Q3 holds for every signed word with all labels in
`{1,3}` and exactly one label `3`, of arbitrary length.

**Proof.**  If the `3` is minus, `J = I(1,c,r)` with `c = #(plus 1)` and
`2r = #(minus 1) + 1`.  If the `3` is plus, the twist `y -> -y` gives
`J = I(1, #(minus 1), (#(plus 1)+1)/2)`.  When the total label parity is
odd, `J = 0` identically.  QED.

The repository covers these words only up to seven factors, and via
Prop. 24A the pattern "minus `{3,1}` plus a fundamental plus suffix".  This
is local acceptance of an infinite family, not theorem-slot closure.
Venus also reduced odd `b >= 3` exactly to a one-variable hypergeometric
integral (packet on file); its sign is open.

### Corollary FM19 (one large symmetric power at r=2; FM-T1-CB, luna_max_mars; accepted)

**Citation.**  De Negri--Sbarra, "Gröbner bases of ideals cogenerated by
Pfaffians", Thm. 2.2: for an anti-diagonal lex order the `6 x 6` Pfaffians
form a Gröbner basis, with leading terms `x_(a1 a6) x_(a2 a5) x_(a3 a4)`.
The standard monomials are therefore exactly the 3-nesting-free loopless
multigraphs in each multidegree.  This settles the fifth-pass model, and
with it `m_11 <= m_20`, unconditionally.

**Statement.**  For `M = Sym^s C^4 tensor (C^4)^(tensor a)`, put `u=(a+s)/2`,
`v=(a-s)/2`, `C_0 = binom(a,u)`, and `L = 5m_00 - 3m_11 + m_20`.

- If `a+s` is odd or `s > a+2`, then `L = 0`.
- If `s = a+2`, then `L = 1`.
- If `4 <= s <= a` and `a+s` is even (interior), then

  ```text
  L = C_0^2 (a+1)(a+2)(s+1)(s+2)(s+3) R_s(v) / ((u+1)^2(u+2)^2(u+3)^2(u+4)(v+1)^2(v+2)),
  R_s(v) = 60v^2 + (-12s^2+12s+180)v + s^4 + 2s^3 - 13s^2 + 10s + 120.
  ```

  Its discriminant is `-48(2s^4+16s^3+22s^2-40s-75) < 0`.
- For `s <= 3`, with `a = s+2v` and `v >= 0` (supplied by FM-CHK3 and
  spot-checked by the main agent at `v = 0`):

  ```text
  s=0: 720(2v+1) binom(2v,v)^2 / ((v+1)^2(v+2)^2(v+3)^2(v+4))
  s=1: 2880(2v+3) binom(2v+1,v+1)^2 / ((v+2)^2(v+3)^2(v+4)^2(v+5))
  s=2: 1440(2v+3)(5v^2+13v+10) binom(2v+2,v+2)^2 / ((v+1)^2(v+3)^2(v+4)^2(v+5)^2(v+6))
  s=3: 2880(2v+5)(5v^2+9v+14) binom(2v+3,v+3)^2 / ((v+1)^2(v+4)^2(v+5)^2(v+6)^2(v+7))
  ```

  Values: `L(0,0)=5`, `L(1,1)=3`, `L(2,2)=2`, `L(3,3)=2`.  The boundary
  `s = a+2` (`v = -1`, where the rational form is undefined) is handled
  separately: only `m_20` contributes there, with multiplicity 1.

The derivation expands the `C_2` Weyl alternant in binomials.  The main
agent verified the closed form exactly on all `s <= 9`, `a <= 14` (150
pairs) and re-derived the discriminant.

**Corollary FM19.**  Q3 holds for every word with minus labels `{q,1,1,1}`
(`q` arbitrary) and any number of plus labels `1`.  By the twist it also
holds for the corresponding sign-swapped words.  This is a new infinite
family for long suffixes.

**Assigned next (FM-T1-AB1).**  Q3 for every word in which all labels but
one equal `1`, at all levels `r`.  Exact screen: 847 words with `q <= 8`,
`#(1s) <= 20`, none negative.
- **FM-CHK2 (luna_max_mercury).**  Fifth-pass multigraph formulas and the
  consequence `m_11 <= m_20` accepted; Lemma FM17 and its gadget
  interpretation accepted.  The `{1,3}` reduction needed count definitions
  and a parity clause (repaired above).
- **FM-T1-SX (luna_max_saturn).**
  - The four-minus family with a fundamental plus suffix reduces exactly to
    a finite fusion-ballot inequality `FB_4` (subset expansion with ballot
    kernels `K_a(l,m) >= 0`).
  - Proved subcases: `a=0` (inside the proved range) and labels forming two
    equal pairs (pointwise nonnegative).
  - The integrand is not pointwise nonnegative: `q=(1,1,1,3)` at
    `(x,y)=(1,0)`.
  - The twelve endpoint values agree with positivity.  No new family, so no
    criterion movement.
  - Next assignment FM-T1-2L: two large symmetric-power parts, with minus
    labels `{s+1,t+1,1,1}` plus a fundamental plus suffix, by the FM19
    Weyl-alternant method.

### Corollary FM20 (three and five 3s; FM-T2b2, luna_max_venus; accepted)

For all integers `b,m,r >= 0`,

```text
I(b,2m,r)/M(m,r) = L_b(m,r)
  = sum_(i+j+k=b) b!/(i!j!k!) 3^i (-2)^k A_i(m) A_j(r)
                                   / ((m+r+2)_(i+j) (m+r+3)_(i+j)),
A_i(s) = prod_(t<i) (2s+2t+1)(2s+2t+3).
```

- For `b=3`, `L_3 = 8 R_3/Q_3`, with `R_3` of degree 6 and
  `Q_3 = (m+r+2)(m+r+3)^2(m+r+4)^2(m+r+5)`.
- For `b=5`, `L_5 = 32 R_5/Q_5`, with `R_5` of degree 10.

Both numerators have positivity certificates on the integer quadrant:

- `R_3(r+t,r)` has all coefficients `>= 64`.  The coefficients of
  `R_3(m,m+z+1)` in powers of `m` are positive-definite quadratics plus
  nonnegative terms for `z >= 0`.
- `R_5(r+z,r)`, `R_5(m,m+t)` for `t = 0,...,4`, and `R_5(m,m+y+5)` have all
  coefficients positive.

The main agent verified symbolically: formula (1) at `b=3` exactly; `R_5`
re-derived from (1) independently of the packet's display; all certificate
coefficient minima (64 and 1024 for `R_5`); and the spot values
`I(5,0,0) = 292`, `I(5,0,1) = 44`, `I(5,4,3) = 1520`,
`I(5,8,6) = 249966`.

**Corollary FM20.**  Q3 holds for every `{1,3}`-word in which all `3`s have
the same sign and the number of `3`s lies in `{0,...,6}` or is even (even
counts are pointwise trivial).  The only open case in this family is an odd
number `>= 7` of `3`s.  For that case, `L_b = E_(m,r)[H_3^b]` with
`H_3 in [-2,10]`, and the main-agent scan shows positive minima (0.74 at
`b=7`, growing fast).  It is assigned as FM-T2b3.

### Kill K14 (two arbitrary plus labels with a fundamental suffix; FM-T2-SX, luna_max_neptune)

Write `T(p,q;a) = <A_1^a B, A_p A_q B>`, where `B = E_(1,0)` and `A_r` is the
additive compound on wedges.  Prop. 24A of the search note keeps positive
wedge vectors positive.  It does not apply here, because `A_p A_q B` has
signed coordinates.  At `(p,q) = (2,4)`:

```text
A_2 A_4 B = E70 - E61 + E52 + 2E50 - E43 - 2E41 + E32 + 2E30 - E21 + E10.
```

Grouping the pairing with `A_1^4 B` by fixed `u+v` gives layer contributions
`3, 7, -2, 0`, so layers cannot be positive termwise.  The exact screen
`p,q <= 4, a <= 6` passes (the first case beyond seven factors,
`(2,4,4)`, has margin 8).  The family remains open.

### Reformulation (main agent): the tilted law is a two-particle Jacobi ensemble

Put `x = 2cos(theta)`, `y = 2cos(phi)`, `s = cos(theta+phi)` and
`t = cos(theta-phi)`, ordered so that `s < t`.  Then

```text
u^2 = 4(1+s)(1+t),   v^2 = 4(1-s)(1-t),
S_2 = 2 + 4st,   H_3 = 1 + (2s+1)(2t+1),   S_3 = u (1 + (2s-1)(2t-1)),
```

and `u^(2m) v^(2r) dsigma(x) dsigma(y)` becomes, up to a constant, the
two-particle Jacobi unitary ensemble with weight
`(1+s)^(m-1/2) (1-s)^(r-1/2)` and Vandermonde `(s-t)^2`.  At `m = 0`,
`r = 1` this is Haar measure on `SO(5)`: the map `(theta,phi) ->
(theta+phi, theta-phi)` is the isogeny `Sp(4) -> SO(5)` on maximal tori.

Matching the FM18 moments gives an independent-Beta form:

```text
u^2 = 16 R C B,   v^2 = 16 R (1-C)(1-B),
R ~ Beta(m+r+1, 1),   C ~ Beta(m+1/2, r+1/2),   B ~ Beta(m+3/2, r+3/2),
```

with `R`, `C`, `B` independent.  All variables are bounded, so the moments
determine the law.

Consequence: Q3 for all words with labels in `{1,2}` is equivalent to
`E_(m,r)[S_2^b] >= 0` for all `m, r, b`.  This family is already a theorem:
Cor. 24B17B and Prop. 24B3 of the search note.  It therefore serves as a known
input for the `{1,2,3}` target below.  With labels in `{1,2,3}` it is
equivalent to `E_(m,r)[S_2^a G^g H_3^d] >= 0` for `g <= 2m`, `d <= 2r`,
where `G = S_3/u`.

Scan: no negative value in 215 exponent triples up to 5 with `m, r <= 10`.
For `S_2^b` the minima are 1, 4, 20.9, 89 at `b = 1, 3, 5, 7`, all on the
diagonal `m = r`.  The `{1,2,3}` family is assigned as FM-L123.

### Corollary FM21 (one plus 3 beside one plus 1; FM-T2-123, luna_max_jupiter; accepted)

Put `tau_(m,b) = <Shat_1^m Shat_2^b, 1>_(Sp4)`.  From
`S_3 = (3 S_1 S_2 - S_1^3 + 2 S_1)/2`, the `Sp(4)` pairing
`<Shat_1^a Shat_2^b Shat_3, 1>_(Sp4) = Phi_00 + Phi_20 - Phi_11` is given
below.  The Q3 word value is `J = ell_1 = 2x` this pairing (normalisation
repaired after FM-CHK4).

```text
(3 tau_(a+1,b+1) - tau_(a+3,b) + 2 tau_(a+1,b))/2.
```

Since `(x+y)^2 = S_2 + 2 + 2xy` and `xy (x^2-y^2)^2` is odd in `x`,

```text
tau_(4,b) = tau_(2,b+1) + 2 tau_(2,b).
```

So at `a = 1` the pairing equals `tau_(2,b+1)`, and `J = 2 tau_(2,b+1)`.
This is `>= 0` by Cor. 24B17B.
Hence Q3 holds for `[1^-, 1^-, 1^+, 2^+ x b, 3^+]` for every `b`.  The main
agent re-derived the packet's table exactly; the values are
`(1,1,1) -> 2`, `(1,2,1) -> 4`, `(3,1,1) -> 4`, `(1,1,3) -> 19`,
`(2,2,2) -> 15`.  It also checked the identity for `b <= 4`, getting
`1, 2, 4, 10, 27`.  For odd `a >= 3` the packet's frontier is
`tau_(a+3,b) <= 3 tau_(a+1,b+1) + 2 tau_(a+1,b)`, which is equivalent to
`E_(m,1)[S_2^b G] >= 0`.  That case is now part of FM-L123.

A further fact from the main agent, checked numerically: in the Jacobi
coordinates, `S_(2n) = 2 + 4 sum_(j=1)^n T_j(s) T_j(t)`, where the `T_j` are
Chebyshev polynomials of the first kind.

### FM-CHK3 (luna_max_mercury, independent checker)

- **FM18: ACCEPTED.**  The checker gave an independent proof of `M(m,n)`.
  Semicircle integration by parts gives
  `(2m+4n+7)A + (2m+1)B = 16(2m+1)M` and the symmetric equation, with
  determinant `8(m+n+2)(m+n+3)`.  The b=1 form, its positivity, and both sign
  placements of the 3 all check.
- **FM19: ACCEPTED** for general `s`.
  - The alternant expansion into the `C_d` quadratic form, the ratio
    substitution, and the symbolic identity with the factorised form all
    check.
  - The discriminant is negative for `s >= 4`.
  - Repairs applied: the interior domain is `4 <= s <= a`; the `s <= 3`
    table is now included; the `s = a+2` boundary is handled separately.
- **FM16 (Grassmannian identity): ACCEPTED.**
  - Atiyah--Bott applicability, the character `sum_I t_I^(-n)`, the `6/24`
    normalisation, and the `u`-constant-term step all check.
  - Repair applied: a full bibliographic pointer.  The checker cited paper I
    (1967); the main agent gave Thm. 4.12 of paper II (1968), the
    holomorphic Lefschetz theorem.  The bibliographic detail is still to be
    confirmed.
- Small checks: `ell_2(1) = 10`, `ell_2(C^4) = 0`, `ell_2(Sym^2) = 2`,
  `ell_2((C^4)^(x)2) = 6`.

### Corollary FM22 (same-sign 3s when the 1s on one side dominate; FM-T2b3, luna_max_venus; accepted)

**Claim.**  For every odd `b` and all integers `m >= r >= 0`,
`I(b,2m,r) >= 0`.

**Proof** (packet; the main agent checked every step and the integral
formula numerically).

1. *Change of variables.*  Put `x = (u+v)/2` and `y = (u-v)/2`, then `v = tu`
   on the first quadrant, then `z = u^2(1+t)^2/16`.  This gives

   ```text
   I(b,2m,r) = 16^(N+1)/pi^2 int_0^inf t^(2r)/(1+t)^(2N+2) J_(N,b)(alpha(t),kappa(t)) dt,
   J_(N,b)(alpha,kappa) = int_0^1 z^N sqrt((1-z)(1-kappa z)) (alpha z - 2)^b dz,
   alpha(t) = 4(3+t^2)/(1+t)^2,   kappa(t) = ((1-t)/(1+t))^2,   N = m+r.
   ```

2. *Stochastic dominance.*  The normalised `z`-law has likelihood ratio
   `sqrt((1-kappa z)/(1-z))`, which is increasing, against `Beta(N+1,2)`.
   `Beta(N+1,2)` in turn dominates `Beta(2,2)` when `N >= 1`.  For odd `b`,
   `(alpha z - 2)^b` is increasing in both `z` and `alpha`.  So
   `F_(N,kappa)(alpha) >= F_*(alpha)`, where `F_*` is the `Beta(2,2)`
   expectation.

3. *The `Beta(2,2)` bounds.*  Put `X = 2z-1`, which is symmetric about 0.
   Then for all `d >= 0`:

   ```text
   F_*(4+2d) = sum_j binom(b,2j) d^(b-2j) (2+d)^(2j) E X^(2j) >= 0,
   F_*(4+2d) + F_*(4-2d) >= 0.
   ```

4. *Pairing `t = s` with `t = 1/s`.*  The two points share `kappa`, and
   - `alpha_h = alpha(s) >= 4`,
   - `alpha_l = alpha(1/s)` lies in `[3,4]`, with `alpha_l >= 3` because
     `(3s-1)^2 >= 0`,
   - `alpha_h + alpha_l = 16(1+s^2)/(1+s)^2 >= 8`.

   The outer weights are `s^(2r)` on the `alpha_h` point and `s^(2m)` on the
   `alpha_l` point.  When `m >= r`, `s^(2r) >= s^(2m)`, and then

   ```text
   s^(2m)(F_h + F_l) + (s^(2r) - s^(2m)) F_h >= 0.
   ```

5. *The case `N = 0`* is a genuine module.

**Check.**  The formula in step 1 matches the exact moment values
`I(7,0,1) = 1222`, `I(7,2,1) = 12576`, `I(7,4,3) = 24400`, `I(3,2,2) = 18`
and `I(5,0,2) = 30` to 1e-9.  The exact scan over odd `b <= 41` and
`r <= m <= 12` found no negative value.

**Word form.**  Take a `{1,3}`-word whose 3s are all minus.  Q3 holds whenever
the number of plus 1s is at least the number of minus labels.  By the twist,
the same holds with plus and minus exchanged.  Together with FM20, the only
open case of the same-sign `{1,3}` family is odd `b >= 7` with `m < r`.

### Lemma FM23 (one arbitrary label, coefficient form; FM-T1-AB1, luna_max_mars)

Define

```text
g_n = [z^(L+n)] (1+z)^a (1-z)^(2r),   L = (a+2r)/2,
h_k = [z^k] (1+z)^a (1-z)^(2r-1),
Delta_m = g_m^2 - g_(m-1) g_(m+1).
```

Then:

```text
W_r(Sym^s C^4 (x) (C^4)^(x)a) = h_(L+s/2)(g_(s/2) - 2g_(s/2+1) + g_(s/2+2)) + g_(s/2+1)(g_(s/2) - g_(s/2+1)),
W_r(Shat_(2m) (x) (C^4)^(x)a) = Delta_m - Delta_(m+1).
```

The derivation writes the torus in the coordinates `X = zw` and `Y = z/w`,
which are the Jacobi coordinates.  There `S_1^a D_1^(2r)` factorises as
`G(X) G(Y)`, and Weyl integration yields 2x2 Toeplitz minors.

The main agent verified both identities exactly: the second for
`r = 1,2,3`, even `a <= 8` and `m <= 4`; the first for even `s <= 6`.

The `g_n` are Krawtchouk values, `g_n = K_(L+n)(2r; a+2r)`.  By Krawtchouk
self-duality, `g_n` is `binom(N, N/2+n)` times an even polynomial of degree
`2r` in `n`.

**`Delta_m >= 0`** (argument made explicit after FM-CHK5).  The `g_n` are the
coefficients `c_k` of `prod_i (1 + a_i z)` with real `a_i` in `{1,-1}`.
Newton's inequalities `E_k^2 >= E_(k-1) E_(k+1)`, where
`E_k = c_k/binom(n,k)`, hold for real `a_i` of any sign; only Maclaurin's
inequalities need positivity.

- If `c_(k-1) c_(k+1) <= 0`, then `c_k^2 >= c_(k-1) c_(k+1)` trivially.
- Otherwise `E_(k-1) E_(k+1) > 0`, and log-concavity of binomials gives
  `c_k^2 >= binom(n,k)^2 E_(k-1) E_(k+1)
         >= binom(n,k-1) binom(n,k+1) E_(k-1) E_(k+1) = c_(k-1) c_(k+1)`.

Check: 4,515 exact instances (`a < 30`, `r <= 7`), no violation.  The
checker's example `(1,0,-2,0,1)` satisfies it: `0 >= -2` and `4 >= 0`.

Parity clauses: the Sym formula needs `a+s` even (otherwise `W_r = 0`); the
`Shat_(2m)` formula needs `a` even and `m >= 1` (odd `a` gives 0); both need
`r >= 1`.  The families need the stronger monotonicity
`Delta_m >= Delta_(m+1)`, and the analogous inequality for Sym.  Exact
screens at `r = 3` found no negative value for `s, a <= 60` (1,861 pairs)
or for `q, a <= 60` (930 pairs).  Status: conditional.

### Two plus labels (FM-T2-PQ, luna_max_jupiter; correct, partial)

- **Decomposition** (verified by the main agent):

  ```text
  T(p,q;a) = sum_(j=|p-q|, step 2)^(p+q) T_1(j;a) + M(p,q;a),
  M(p,q;a) = int U_p(x) U_q(y) (x+y)^a (x-y)^2 dsigma dsigma.
  ```

  `M` has an exact finite ballot-number formula and vanishes for
  `p+q > a+2`.
- **Proved:** `p = q` for every `a` (the integrand is a square times an even
  power); `a = 0` for all `p, q`.
- **Conditional:** `p+q > a+2`, on the `r = 1` single-label statement
  `T_1(j;a) >= 0`.
- **Open:** `p != q` with `p+q <= a+2`, where `M < 0` occurs, e.g.
  `(2,4,4)` has `M = -1`.

### Two large symmetric powers at r=2 (FM-T1-2L, luna_max_saturn; correct, conditional)

Take `L(s,t,a) = W_2(Sym^s (x) Sym^t (x) (C^4)^(x)a)`.  The packet writes `L`
exactly as a signed sum of chamber-walk counts `K_a`, which the reflection
principle gives in closed form.  It proves the following cases:

| case | value of `L` |
|---|---|
| `s+t+a` odd | 0 |
| `|s-t| > a+2` | 0 |
| `|s-t| = a+2` | 1 |
| `t = 0` | FM19 |
| `s = t` | `>= 0` (a square) |
| `a = 0` | 5, 3, 1 or 0 |

Exact formula (packet FM-T1-2L, formula (1); inserted after FM-CHK5).  Let
`K_a(x,y)` count length-`a` walks with unit steps `+-e_1, +-e_2` from `x` to
`y`, staying in the open chamber `x_1 > x_2 > 0` (shifted coordinates
`lambda + rho`, `rho = (2,1)`).  By reflection,

```text
K_a(x,y) = sum_(w in W(C_2)) det(w) P_a(x - w y),
P_a(d_1,d_2) = sum_k binom(a,k) b_(a-k)(d_1) b_k(d_2),
b_n(d) = binom(n,(n+d)/2) if |d| <= n and n+d is even, else 0.
```

With `x_ij = (s+t-j-i+2, j-i+1)`, summing over `0 <= i <= j <= t`
(`Sym^s (x) Sym^t = sum_(i<=j<=t) V_(s+t-j-i, j-i)`, `s >= t`),

```text
L(s,t,a) = sum_(j=0)^t sum_(i=0)^j [5K_a(x_ij,(2,1)) - 3K_a(x_ij,(3,2)) + K_a(x_ij,(4,1))].
```

FM26 uses the equivalent 24-term Weyl-binomial form of each summand
`Lambda(p,q;a)`, which the main agent re-derived and checked independently.

The first interior case is `L(2,1,1) = 2`.  The constituent `V_(2,1)`
contributes `-2`, so the terms are not individually positive.  The main agent
re-derived all 20 screen values exactly.  The interior `s > t >= 1`,
`a >= 1`, `s-t < a+2` is open and is assigned as FM-T1-2L2, using
closed-form Weyl alternants and telescoping.

### FM-CHK4 (luna_max_mercury, independent checker)

- **Jacobi reformulation: ACCEPTED.**
  - Checked: the fundamental domain `0 <= phi <= theta`,
    `theta + phi <= pi`; the pushforward density
    `(2/pi^2)(s-t)^2/sqrt((1-s^2)(1-t^2))`; the claim that the four-point
    fibres flip the signs of `u` and `v` independently, so odd powers
    average to zero; and the SO(5) Haar identification via the `B_2`
    root product.
  - Boundary check at `theta = phi = pi/3`: `u^2 = 4`, `v^2 = 0`, `S_2 = 0`,
    `H_3 = 1`, `S_3 = -2`.
- **Independent-Beta form: ACCEPTED.**  The moments match, and moment
  determinacy holds on `[0,16]^2`.  Small cases: `E_(0,1)[u^2] = 1`,
  `E_(0,1)[v^2] = 5`, `E_(0,1)[u^2 v^2] = 3`.
- **Chebyshev identities: ACCEPTED.**  `S_(2n) = 2 + 4 sum_(j=1)^n T_j(s) T_j(t)`,
  and for odd labels `S_(2n+1)/u = sum_(j=0)^n P_j(s) P_j(t)`, where
  `P_0 = 1`, `P_1 = 2z-1`, `P_(j+1) = 2z P_j - P_(j-1)`, so that
  `P_j(cos 2g) = cos((2j+1)g)/cos g`.
- **FM20: ACCEPTED.**  The general-`b` formula is re-derived, and the
  `R_3` and `R_5` certificates are re-checked independently.  `R_3` is
  displayed in full in the packet.  The limits are `L_3 -> 8` and
  `L_5 -> 32` as `r -> inf`.
- **FM21: DEFECT, now repaired.**  The displayed formula is the `Sp(4)`
  pairing, which is half the Q3 word value `J`.  The recursion and the
  positivity are unaffected.
- **Atiyah--Bott citation: ACCEPTED.**  Part II (Ann. of Math. 88 (1968),
  451--491), Thm. 4.12, is the holomorphic Lefschetz formula (checker
  confidence 0.98).  The "to be confirmed" flag in FM-CHK3 is cleared.

### Corollary FM24 (one arbitrary label at r=1; FM-T2-PQ2, luna_max_jupiter; accepted)

Set `N = a+2`, `L = N/2` and `t = m(m+1)`.  Let `m` be an integer for even
`a`, a half-integer for odd `a` (the odd labels), with `0 <= 2m <= N`.  Then

```text
T_1(2m;a) = <Shat_(2m) (C^4)^(x)a, 1>_(Sp4)
          = binom(N,L+m)^2 (2m+1) (3N^2 + (6-8t)N + 16t(t-1))
            / (N(N-1)(L+m+1)^2(L+m+2)(L-m+1)).
```

This follows from FM23 at `r = 1` with
`g_n = binom(N,L+n)(4n^2-N)/(N(N-1))`, the Krawtchouk `K_2` form.

The last factor is positive in every case:

- at `m = 0` it is `3N^2 + 6N`;
- at `m = 1/2` it is `3N^2 - 3`;
- for `m >= 1` (so `t >= 2`) its discriminant `-128t^2 + 96t + 36` is
  negative.

**Checks.**  The main agent verified the closed form symbolically (`sympy`,
identity in `N` and `m`), and against direct integration for all 460
parity-compatible `(j,a)` with `a < 40`.

**Corollary FM24.**  Q3 holds for every word `[1^-, 1^-, j^+, 1^+ x a]`, and
by the twist for its sign-swapped counterpart.  With the two-label
decomposition this gives, unconditionally:

- `T(p,q;a) >= 0` whenever `p+q > a+2`;
- `T(1,q;a) = T_1(q;a+1) >= 0`;
- the cases `p = q` and `a = 0`.

**Telescoping (main agent, verified exactly).**  Since
`T_1(j;a) = Delta_(j/2) - Delta_(j/2+1)`,

```text
T(p,q;a) = Delta_(|p-q|/2) - Delta_((p+q)/2+1) + M(p,q;a).
```

Here `M(p,q;a) = <V_p x V_q, S_1^a D_1^2>_K`.  By the quadrant reflection
principle it is a short signed sum of products of two binomials.  Both
identities were checked exactly: the telescoping for `a < 16`,
`p, q <= 8`, and the reflection form of `M` for `a < 12`, `p, q <= 7`.  The
open regime `p != q`, `p+q <= a+2` is therefore a closed-form, three-parameter
binomial inequality (assigned FM-T2-PQ3).  Jupiter's packet adds that the
unrestricted quadratic form in the `g_n` is not positive semidefinite, e.g.
direction `(1,1,1)` at `(2,4,4)`, so a proof must use the specific `g_n`.

### Corollary FM25 (labels {1,2,3}, partial; FM-L123, luna_max_neptune; accepted)

Work in the Jacobi form with `S = S_2 = 2+4st`, `T = s+t`, `N = m+r` and
`Delta = m-r`.  Then `G = S - 2T` and `H_3 = S + 2T`.  Put
`tau_j = E[S^j]`, which is `>= 0` by Prop. 24B3 / Cor. 24B17B, and
`rho_j = E[T S^j]`.  Integration by parts against `h = (1-z^2) w` gives

```text
(N+j+3) rho_j = 2 Delta tau_j + 6 j rho_(j-1),
(N+j+4) E[T^2 S^j] = 2 Delta rho_j + tau_j + tau_(j+1)/2 + 6 j E[T^2 S^(j-1)].
```

By induction `rho_j` has the sign of `Delta`.  Hence:

- `E[S^a G] = tau_(a+1) - 2 rho_a >= 0` for all `a` when `m <= r`;
- `E[S^a H_3] = tau_(a+1) + 2 rho_a >= 0` for all `a` when `m >= r`.

In addition, exact closed forms with all-positive polynomial numerators
settle:

- `a = 1, 2`, all `m, r`, for one 3 (FM18 already covers `a = 0`);
- the mixed moments `E[G H_3]` and `E[S G H_3]` for `m, r >= 1`.

**Checks.**  The main agent verified both recurrences exactly for
`m, r <= 6`, `j <= 5`, and the four closed forms `E[SG]`, `E[S^2 G]`,
`E[G H_3]` and `E[S G H_3]` exactly on grids up to 8.

**Corollary FM25.**  Q3 holds for every `{1,2,3}`-word with exactly one 3 in
the following cases:

- the 3 is plus and `2m <= 2r`, where `2m = #1^+ + #2^- + 1` and
  `2r = #1^- + #2^-`;
- the 3 is minus and `2m >= 2r`, with the analogous counts;
- there are at most two plus 2s, with any counts.

It also holds for every `{1,2,3}`-word with one plus 3, one minus 3, at most
one plus 2, and `m, r >= 1`.

**First open case:** `a = 3`, `m > r`, one plus 3, i.e.
`tau_4 >= 2 rho_3`.

### Corollary FM26 (two symmetric powers, smaller one at most 2; FM-T1-2L2, luna_max_saturn; accepted)

Write `L(s,t,a) = W_2(Sym^s (x) Sym^t (x) (C^4)^(x)a)`.  The interior of the
two slices has the closed forms

```text
L(s,1,s-1+2v) = C_1^2 (s+1)(s+2)(s+3)(s+2v)^2(s+2v+1)(s+2v+2) R_s(v)
                / ((s+v)^2(v+1)^2(v+2)(s+v+1)^2(s+v+2)^2(s+v+3)^2(s+v+4)),
L(s,2,s-2+2v) = C_2^2 (s+1)(s+2)(s+3)(s+2v)(s+2v-1) Pi_2(s,v)
   / ((s+v)^2(v+1)^2(v+2)(s+v-1)^2(s+v+1)^2(s+v+2)^2(s+v+3)^2(s+v+4)),
```

where

- `C_1 = binom(s-1+2v, s+v-1)` and `C_2 = binom(s-2+2v, s+v-2)`;
- `Pi_2` is given explicitly by (with `n = s-3`)

  ```text
  Pi_2 = (n+2)(n+3)(n+4)(n+5)(n+6)(n+7)(n^2-n+4)
       + 2(3n^7+61n^6+479n^5+1965n^4+5646n^3+15026n^2+28660n+22992) v
       + 2(8n^6+120n^5+675n^4+2426n^3+9037n^2+23566n+23880) v^2
       + 4(5n^5+35n^4+51n^3+757n^2+5296n+8562) v^3
       + v^4 (A(n) + B(n) v + 600 v^2),
  A(n) = 2(5n^4-50n^3-149n^2+2810n+8076),   B(n) = -120(n^2-5n-37),
  B^2 - 2400A = -4800(2n^4-20n^3-2n^2+1700n+3969) < 0;
  ```

  here `A(n) > 0` for all `n >= 0` (checked piecewise in the packet);
- `R_s` is the same quadratic as in FM19;
- `Pi_2 > 0` has an explicit certificate: its coefficients are positive, and
  the leftover quadratic in `v` has negative discriminant.

The derivation is a 24-term Weyl-binomial formula for
`Lambda(p,q;a) = W_2(V_(p,q) (x) (C^4)^a)`, summed over the constituents of
`Sym^s (x) Sym^t`.  The main agent re-derived `Lambda` independently and
checked both closed forms exactly for `s <= 11`, `v <= 7`.

**Corollary FM26.**  Q3 holds for minus labels `{s+1, t+1, 1, 1}` with any
number of plus 1s whenever `min(s,t) <= 2`.  The case `s > t >= 3` is open.
The recurrence of `R_s` across `t = 0, 1` suggests a general-`t` factorisation.

### FM-T2b4 (luna_max_venus): the remaining `{1,3}` case as a Hankel sum

Let `X ~ Beta(m+1/2, r+1/2)`, `xi = 4X - 1` and `mu_j = E xi^j`.  Then

```text
I(b,2m,r)/M(m,r) = sum_k binom(b,k) D_k / D_0,   D_k = mu_k mu_(k+2) - mu_(k+1)^2
```

(Andreief).  The packet also gives an exact obstruction: the `Beta(2,2)`
comparison used for FM22 cannot handle `m < r`, since the paired bound is
negative at `b = 7`, `(m,r) = (0,1)`, `s = 1/19999`.  The exact law repairs
this to first order.  Still open.

## Mechanism investigation (main agent, 2026-09-27)

Goal: a structural reason for Q3 that covers every word, as opposed to one
closed form per family.

**M0 (Ginibre rotation).**  Put `alpha = theta+phi`, `beta = theta-phi` and
`K_p = {-p, -p+2, ..., p}`.  Then

```text
S_p = 2 sum_(k in K_p) cos(k alpha/2) cos(k beta/2),
D_p = -2 sum_(k in K_p) sin(k alpha/2) sin(k beta/2).
```

So a word with `n` factors, `d` of them minus, equals
`2^n (-1)^d sum_kappa F_kappa(alpha) F_kappa(beta)`.  For even `d` this is a
positive-semidefinite kernel.  For odd `d` the full integral vanishes by the
swap.  (Scalar repaired after FM-CHK6.)  Here
`F_kappa = prod_i cos(k_i alpha/2)` or `sin(k_i alpha/2)`.  The Haar weight
is `sin^2 theta sin^2 phi = (cos alpha - cos beta)^2/4`, the SO(4) Weyl
density.  Hence

```text
J = 2^n sum_kappa Q(F_kappa)  (d even),   Q(f) = f_0^2 + f_0 f_2 - 2 f_1^2
                                 = 2(<F><F c^2> - <F c>^2),   c = cos alpha.
```

- For U(1) the weight is 1, so `J = sum (int F_kappa)^2`.  This is Ginibre's
  proof for plane rotators.
- For SU(2) each term is a 2x2 Gram determinant of a signed weight, and `Q`
  has Lorentz signature `(1,2)`.
- Terms are not individually positive.  For `S_2`, the magnitude-2 term
  gives `-2` and the magnitude-0 term gives `+2`.
- The grouping into full characters (all of `K_p`) is essential: a single
  term `cos alpha (x) cos beta` gives `J < 0`.
- Combinatorial form:
  `J = sum_(magnitudes) c [A(0)^2 + A(0)A(4) - 2A(2)^2]`, where `A` is the
  signed count of `sum_i sigma_i m_i`.

**M1 (moment form and boundary positivity).**

- *Coordinates.*  Put `a = u^2/16` and `b = v^2/16`.  These lie in
  `Omega = {sqrt(a) + sqrt(b) <= 1}`, and `(a,b) = (RCB, R(1-C)(1-B))`.
- *Equivalence.*  Fix the non-1 part `P` of a word.  Q3 for every
  1-suffix is equivalent to `int a^m b^r P dnu_0 >= 0` for all `m, r`, with
  half-integers allowed.
- *The boundary.*  The outer boundary `sqrt(a) + sqrt(b) = 1` is the set
  `{x = +-2 or y = +-2}`, where one SU(2) factor is central.  At `x = 2`
  every generator is nonnegative: `S_p = p+1 + U_p(y)` and
  `D_q = q+1 - U_q(y)`.
- *Asymptotic positivity.*  For large `(m,r)` the weight `a^m b^r`
  concentrates on this boundary, which explains positivity asymptotically.
  The density vanishes there to second order (Vandermonde), which matches
  the asymptotic tightness of T1.
- *What remains.*  Finite `(m,r)` is the real difficulty.

**Kill K15 (orthant dominance).**  The candidate: for every upper set `U`
in the `(u,v)` quadrant, `int_U P dsigma dsigma >= 0`.  This would give all
moments at once, but it is false.  Normalised minima on a 300-point grid:

| `P` | minimum |
|---|---|
| `D_4` | -0.33 |
| `S_5`, `D_5` | -0.29 |
| `S_4` | -0.22 |
| `(2+,5+)` | -0.10 |
| `S_3`, `D_3` | -0.08 |

All these words satisfy Q3.  The negative mass is oscillatory
(Chebyshev), not stochastically below the positive mass.

**M2 (Rodrigues).**  `U_p dsigma = (-d/dx)^p [c_p (4-x^2)^(p+1/2) dx]`.  For
single labels this gives

```text
J = 2 sum_k binom(p,k) (c_+)_k (c_-)_(p-k) I_k.
```

The terms with `p-k` even are nonnegative.  The terms with `p-k` odd can be
negative (52 cases with `p <= 6`), and are dominated only by both even
neighbours together.  This is again a Turán/log-concavity condition, the
same one that closes FM23 and FM24.

**Mechanisms in the accepted families.**

| mechanism | families |
|---|---|
| explicit closed forms, positivity by discriminant | FM18, FM19, FM20, FM24, FM26 |
| stochastic dominance with twist pairing `t <-> 1/t` | FM22 |
| Jacobi integration by parts with sign propagation | FM25 |
| Toeplitz minors `Delta_m - Delta_(m+1)` (Turán monotonicity) | FM23, FM24 |

The common thread: positivity always comes from a monotonicity (Turán or
stochastic) statement in the radial/Beta variables, after the Ginibre
rotation.  Delegated as FM-MECH-A (certificate classes closed under
multiplication by generators) and FM-MECH-B (positive recursions/semigroups
on words).

**Literature check (main agent, web search 2026-09-27; abstracts only).**

- *Ginibre (1970).*  The rotation mechanism M0 is exactly Ginibre's proof
  for plane rotators.
- *Sylvester, CMP 73 (1980), "The Ginibre inequality".*  The graph
  (cycle-group) Ginibre inequality holds for spin dimensions 1 and 2 and
  fails on some graphs for every dimension `>= 3`.  That result concerns
  O(n) pair interactions `sigma_x . sigma_y`, not the single-site condition
  for SU(2) class functions (zonal functions on `S^3`) studied here, so it
  does not contradict Q3.  Q3 is, however, the `S^3`-zonal analogue, and the
  failure for `n >= 3` in the graph setting shows that dimension-free
  arguments cannot work.  This is consistent with the rank dependence of T1
  and T2.
- *Abdesselam (arXiv 2207.07603).*  Proves GKS2/Ginibre-type inequalities
  for O(N) and CP^(N-1) in an asymptotic regime, via stable determinantal
  (Kirchhoff) polynomials and the Rayleigh property.  This is the nearest
  known "Lorentzian/stable" mechanism.  Compare the Lorentz signature of
  `Q` in M0.
- *Chevyrev--Garban (arXiv 2404.09928).*  The Villain limit covers
  non-abelian groups, but Ginibre is proved only in the abelian case.

No SU(2)-character Ginibre mechanism was found in these sources.

### FM-CHK5 (luna_max_mercury, independent checker)

- **FM22: ACCEPTED in full.**  Checked: density, constant, dominance chain,
  pairing, `N = 0`, and the word translation for both placements of the
  3s (`#S_1 >= #D_3 + #D_1`, or twisted `#D_1 >= #S_3 + #S_1`).
- **FM23: identities ACCEPTED**, with the parity clauses now stated.  The
  checker disputed the Newton step.  The main agent disagrees and has made
  the argument explicit above: Newton holds for real roots of any sign.
- **Two-label decomposition: ACCEPTED.**  Also checked: `M(2,4;3) = 0`,
  `M(2,4;4) = -1`, `T(1,3;0) = 1`.
- **K14: ACCEPTED.**  Vectors and layer values are exact.
- **Saturn section: NEEDS-REPAIR (repaired).**  Formula (1), the kernel,
  the index limits and the boundary conventions are now in the note.  The
  boundary values were re-checked independently.

**Where the cancellation is (main agent, exact values against numerical
`|F|`).**  The data below cover 1,796 signed words with labels `<= 5` and
length `<= 6`.

- *No accidental zeros.*  Among 1,514 words with labels `<= 4` and length
  `<= 6`, every zero is forced by parity or by support (the largest label
  exceeding the sum of the others).  There are no nontrivial equality
  cases.
- *Worst cancellation.*  The smallest ratios `J / int |F| dmu` occur at the
  support boundary, where the largest label equals the sum of the others.
  There `J = 2` comes from the unique top-weight path, yet the integrand is
  large:

  | word | `J` | `int |F|` | ratio |
  |---|---|---|---|
  | `[1^-x5, 5^-]` | 2 | 42.6 | 0.047 |
  | `[1^-x4, 4^+]` | 2 | 16.0 | 0.125 |
  | `[1^-x3, 2^+, 5^-]` | 2 | 12.7 | 0.158 |

- *Consequence.*  Any purely analytic certificate must reproduce these
  highest-weight cancellations exactly.  Combinatorially they are trivial,
  since only the top constituent survives.  Long words with small labels
  are the opposite case: the measure concentrates at the central boundary
  and positivity is analytically easy.  This suggests a hybrid mechanism,
  with highest-weight or support structure near the boundary and analytic
  positivity in the interior.
- *Tightness relative to the all-plus word* is largest for balanced
  `{1,2}` words (ratio 0.2, 0.077, 0.020, 0.0049 at lengths 4, 6, 8, 10).
  Those words are already proved, and pure 1-words are pointwise positive,
  so that ratio is not the difficulty measure.

### Lemma FM27 (one arbitrary label reduces to a single Turán monotonicity; FM-T1-AB2 (mars) + main agent)

Let `c_j = [z^j] (1+z)^a (1-z)^t`, `N = a+t` and `D_k = c_k^2 - c_(k-1) c_(k+1)`.

**Reduction (packet FM-T1-AB2).**  Every one-arbitrary-label family, at every
level `r`, has the form `W_r = D_k - D_(k+1)`:

- `Shat_q (x) (C^4)^(x)a` for either parity of `q`: here `t = 2r` and
  `k = (N+q)/2`;
- `Sym^s C^4 (x) (C^4)^(x)a`: here `t = 2r-1` and `k = (N+s+1)/2`.

Off parity, the value is 0.  So the family's Q3 is exactly the monotonicity
`D_k >= D_(k+1)` on the target indices, all of which satisfy `2k > N`.  The
main agent verified both formulas exactly against direct integration for
`r = 1,2,3`, `a <= 8`, `q <= 7` and `s <= 6` (189 cases, including the zero
values off parity).

**Identity (main agent).**  From the Krawtchouk recurrence
`(k+1) c_(k+1) = (a-t) c_k - (N-k+1) c_(k-1)`:

```text
(N-k+1) D_k - (k+2) D_(k+1) = c_k^2 - c_(k+1)^2.                    (FM27a)
```

Put `omega_j = 1/((N+2) binom(N+1,j+1))`, which is increasing for
`j >= (N-1)/2`.  Summing (FM27a) gives, for `N/2 <= k < N`,

```text
D_k - D_(k+1) = c_k^2/(N-k+1) + c_(k+1)^2 (3k-2N+1)/((N-k+1)(N-k))
              + (binom(N+2,k+1) - binom(N+2,k+2)) sum_(j>=k+2) c_j^2 (omega_j - omega_(j-1)).   (FM27b)
```

**Checks.**  (FM27a) holds in 33,750 exact cases and has a two-line proof
from the recurrence.  (FM27b) holds exactly in every right-half case of the
same grid.

**Consequence (proved).**  For `k >= (2N-1)/3`, every term of (FM27b) is
`>= 0`, so `D_k >= D_(k+1)`.  In word terms, Q3 holds for every word whose
labels are all 1 except one label `L`, at every level `r`, whenever
`L >= (N-2)/3`, with `N` counting the fundamental labels:

- for `Shat_q`, this is `q >= (a+2r-2)/3`;
- for `Sym^s`, it is `s >= (a+2r-6)/3`.

The packet's discriminant condition (*) covers `a ~ t`.  At `r = 1`, FM24
covers everything.

**Open (middle range).**  `N/2 < k < (2N-1)/3` with `|a-t|` large.  An
exact scan (`a <= 160`, `t <= 21`, 147,320 target instances, plus all
right-half `k`) found no violation of `D_k >= D_(k+1)`.

### Theorem FM28 (half-quadrant cone mechanism; FM-MECH-B, luna_max_venus; accepted)

Work in the Jacobi form with `T = s+t`, `S = S_2 = 2+4st`, `N = m+r` and
`Delta = m-r`, and put
`A = (1-s^2) d_s + (1-t^2) d_t`.  Integration by parts against
`h = (1-z^2) w` gives, for every polynomial `f`,

```text
E_(m,r)[ A f + (2 Delta - (N+3) T) f ] = 0,
A T = 1 + S/2 - T^2,   A S = T(6 - S).
```

Hence, with `R_(a,b) = E[T^a S^b]`,

```text
(N+a+b+3) R_(a+1,b) = 2 Delta R_(a,b) + a R_(a-1,b) + (a/2) R_(a-1,b+1) + 6b R_(a+1,b-1).   (FM28)
```

The base case `R_(0,b) = E[S^b] >= 0` is Prop. 24B3 / Cor. 24B17B (the
`{1,2}` theorem).  Induction on `a+b` gives:

- `R_(a,b) >= 0` whenever `m >= r`;
- `(-1)^a R_(a,b) >= 0` whenever `m <= r`.

**Theorem FM28.**

- If `m >= r`, then `E_(m,r)[P(T,S)] >= 0` for every polynomial `P` with
  nonnegative coefficients in the monomials `T^a S^b`.
- If `m <= r`, the same holds with `-T` in place of `T`.

Both cones are closed under multiplication.  This is a genuine semigroup
mechanism, valid on a half-quadrant.

**Which generators lie in the cones** (main agent, exact).

| generator | `(T,S)` form | cone |
|---|---|---|
| `S_1/u`, `H_2/u` | `1` | both |
| `S_2` | `S` | both |
| `H_3` | `S + 2T` | `C+` |
| `S_3/u` | `G = S - 2T` | `C-` |
| `S_4` | `S^2 + S - 8T^2` | neither |
| `H_4/u` | `S - 1` | neither |
| `S_5/u` | `S^2 - 2ST - 4T^2 + 4T - 1` | neither |
| `H_5` | `S^2 + 2ST - 4T^2 - 4T - 1` | neither |

The same holds for every label up to 7: no generator with a label `>= 4` lies
in either cone.

**Corollary (words).**  Take a `{1,2,3}`-word whose 3s all have one sign.

- If all 3s are minus, Q3 holds when `#1^+ >= #1^- + #3^-`.  Here
  `2m = #1^+ + #2^-` and `2r = #1^- + #2^- + #3^-`.
- If all 3s are plus, Q3 holds when `#1^- >= #1^+ + #3^+` (twist).
- The number of 2s of either sign is arbitrary.

This contains FM22 and the half-quadrant part of FM25.

**Checks.**  The main agent verified (FM28) exactly for `m, r <= 5`,
`a, b <= 3`.  The sign pattern holds for `m, r <= 7`, `a <= 5`, `b <= 4`.
`E[S^c H_3^d] >= 0` holds for all `r <= m <= 8`, `c <= 4`, `d <= 7`.

**Limits.**

- The `{1,3}` case with `m < r` (e.g. `E_(0,1)[H_3^7] = 611`) is an
  alternating sum of certified-sign terms.  Its terms are
  `1696, -8498, 18879, -22680, 16940, -7392, 1862, -196`.
- Labels `>= 4` fall outside both cones.

Extending the cone is the next step of FM-MECH-B.

### Conjecture HQ (half-quadrant cone; main agent, exact evidence)

The exact `(T,S)` forms of the generators up to label 8 are:

| label | plus: `S_p` | minus: `H_q` |
|---|---|---|
| 2 | `S` | `u v * 1` |
| 3 | `u * (S - 2T)` | `v * (S + 2T)` |
| 4 | `S^2 + S - 8T^2` | `u v * (S - 1)` |
| 5 | `u * (S^2 - 2ST - 4T^2 + 4T - 1)` | `v * (S^2 + 2ST - 4T^2 - 4T - 1)` |
| 6 | `S^3 + S^2 - 12ST^2 - 2S + 16T^2 - 2` | `u v * (S^2 - S - 4T^2)` |
| 7 | `u * (S^3 - 2S^2T - 8ST^2 + 4ST - 2S + 8T^3 + 8T^2 - 4T)` | `v * (S^3 + 2S^2T - 8ST^2 - 4ST - 2S - 8T^3 + 8T^2 + 4T)` |
| 8 | `S^4 + S^3 - 16S^2T^2 - 3S^2 + 20ST^2 - 2S + 32T^4 - 16T^2 + 2` | `u v * (S^3 - S^2 - 8ST^2 - S + 12T^2)` |

(The `u`, `v` factors of each minus label are included.)

**Evidence** (exact moment formula, validity constraints on `m, r` from the
atoms' `u`, `v` factors):

- For `m >= r`, every generator except the plus-odd labels `S_3, S_5, S_7`
  has `E_(m,r)[Phi T^a S^b] >= 0` for `a, b <= 3` and `m <= 8`.
- The plus-odd labels fail (e.g. `S_3` gives -9.75 at `(m,r,a,b) = (1,1,3,3)`).
  They belong to the twisted cone.
- Products of two and three admissible generators, times `T^a S^b`
  (`a, b <= 1`), `r <= m <= 6`: 28,392 cases, no negative value.

**Conjecture HQ** (convention clarified after FM-MECH-A2).  Take a product

```text
X = T^a S^b prod_i Phi_i,
```

where each `Phi_i` is the `(T,S)` quotient polynomial of a plus-even
generator `S_(2k)` or of a minus generator `H_q`.  These are exactly the
entries in the table above, with the `u`, `v` prefactors removed.  Let
`(m,r)` be the TOTAL exponents of the word, so that `2m` and `2r` include
the `u`, `v` prefactors of the generators.  The conjecture: if `m >= r`,
then `E_(m,r)[X] >= 0`.

The prefactors are counted in `(m,r)` and are **not** multiplied into `X`
again.  Multiplying them in (the reading used in FM-MECH-A2) shifts `r`
upward and leaves the half-quadrant.  For example, `E_(2,2)[T Phi_5^2] = +1`
but `E_(2,3)[T Phi_5^2] = -3/20`.
By the twist, the mirror statement holds for `m <= r`, using `-T` and the
plus-odd labels while excluding the minus-odd labels `>= 3`.

**Word form.**  Q3 would hold for:

- every word with no plus-odd label `>= 3` and `#1^+ >= #1^- + #(odd minus labels >= 3)`;
- every word with no minus-odd label `>= 3` and `#1^- >= #1^+ + #(odd plus labels >= 3)`.

This is the first candidate mechanism that covers unbounded labels.
FM28 is its restriction to labels `<= 3`.  Assigned as FM-MECH-B2.

### Corollary FM29 (labels {1,2,3}, extended; FM-L123b, luna_max_neptune; accepted)

**Sign rule.**  This is neptune's independent route to the FM28 sign
pattern.  Put `mu_(j,k) = E[T^k S^j]`.  Then `mu_(j,k) >= 0` for even `k`,
and has the sign of `Delta` for odd `k`.

**Companion recurrence.**  Integrate by parts with `F = (s-t) f(st)` and
symmetrise:

```text
(N+j+2) tau_(j+1) = 2 Delta rho_j + 2N tau_j + 8j upsilon_(j-1) + 4j tau_(j-1),
(N+a+2) E[S^a G]  = 2N tau_a + 4a tau_(a-1) + 8a upsilon_(a-1) - 2(2r+a+2) rho_a,
```

where `rho_j = E[T S^j]` and `upsilon_j = E[T^2 S^j]`.

**Proved** (in addition to FM25 and FM28):

- *One 3, `a <= 3` two's, full quadrant.*  The new case is
  `E[S^3 G] = 16 P(e,r)/((N+2)(N+3)^2(N+4)^2(N+5)^2(N+6))` with `e = m-r >= 0`
  and `P` having all coefficients positive.
- *One plus 3 on the hard half, odd `a >= 2r+2`.*  This comes from
  Cauchy--Schwarz (`rho_a^2 <= upsilon_(a-1) tau_(a+1)`) together with the
  second recurrence.
- *One plus 3 and one minus 3* (`F_a = E[S^a G H_3]`):
  - `a <= 3` for all `m, r >= 1`, with explicit all-positive `Q_2`, `Q_3`;
  - all `a >= N+1`;
  - even `a` with `N <= 2a`.

**Checks.**  The main agent verified exactly for `m, r <= 6`: the companion
recurrence, the `E[S^a G]` identity, `E[S^3 G]`, `F_2` and `F_3`.

**First open formulas.**

- `E[S^4 G] >= 0` for `m > r` (value 28 at `(m,r) = (1,0)`);
- `E[S^4 G H_3] >= 0` for `m+r >= 9`.

These words lie outside both cones of Conjecture HQ (wrong half-quadrant, or
3s of both signs).  Assigned as FM-L123c.

### Kills K16–K19 (FM-MECH-A, luna_max_saturn; certificate route)

- **K16.**  The derivative cone of positive, separately even measures
  supported on `Omega` is not closed under `S_2`.  It contains `delta_0`,
  but `S_2 delta_0 = -2 delta_0` has negative mass.
  - Positive, however: every pure-1 word has an explicit certificate of order
    `<= 2` (tail integrals), and so does `D_2 = uv`.
- **K17.**  The unpaired x/y-Rodrigues construction followed by rotation
  fails at `D_2`.  The coefficient of `d_u^2` is
  `C rho_1(1) rho_1(0)(9-16) < 0` at `(x,y) = (1,0)`.
- **K18.**  Grouping `J = sum Q(F_kappa)` by frequency-magnitude orbits is not
  termwise positive at `S_2`: the orbit contributions are `+2` and `-2`.
  So an orbit is not a Gram determinant of a positive measure.
- **K19.**  The FM23 coefficient sequence is not PF, even after a
  checkerboard twist.  At `a = 2, r = 2`:
  `g = (1,-2,-1,4,-1,-2,1)`, with twist `(1,2,-1,-4,-1,2,1)`.  Its Turán
  minors `9, 5, 1` still decrease.

Conclusion: the live mechanism is the recursion cone (FM28 and Conjecture
HQ), not certificate cones.

### Corollary FM30 (two plus labels: closed form and boundary; FM-T2-PQ3, luna_max_jupiter; accepted)

Take `p > q >= 2`, `a+p+q` even and `p+q <= a+2`.  Put `N = a+2`,
`d = (p-q)/2`, `e = (p+q)/2` and `beta_x = binom(N, N/2+x)`.  Then

```text
M(p,q;a) = beta_(e-1) beta_(d-1) 16 (N-2e+2)(p+1)(q+1) P_M
           / (N(N-1)(N+2d)(N+2e)(N+2d+2)(N+2e+2)(N+2e+4)),
P_M = 3N^2 - 4N d^2 - 4N e^2 - 8N e + 6N + 16 d^2 e^2 + 32 d^2 e + 8 d^2 - 8 e^2 - 16 e,
T(p,q;a) = Delta_d - Delta_(e+1) + M.
```

The main agent verified the closed form for `M` exactly against the ballot
formula: 381 cases, `a < 24`, `p <= 13`.

**Boundary `p+q = a+2` (proved).**  With `s = p-q` and `beta = binom(N,p)`,

```text
T = beta/(N^2(N-1)) [beta R/((p+1)(q+1)) + N(s^2-N)],   R = s^4 + 14pq - p^2 - q^2.
```

This is positive because `R - 2(p+1)(q+1) > 0` and `beta >= N(N-1)/2`.

**Corollary FM30.**  For `[1^-, 1^-, p^+, q^+, 1^+ x a]`, Q3 now holds in all
of the following cases:

- `p+q >= a+2`;
- `p = q`;
- `p = 1` or `q = 1`;
- `a = 0`.

**Open (strict interior, `h = (a+2-p-q)/2 >= 1`).**  The remaining case is
the exact one-variable inequality

```text
delta_d + gamma rho - delta_(e+1) rho^2 >= 0,
rho = prod_(i=0)^q (h+i)/(p+q+h+1-i),
```

with `delta` and `gamma` explicit (see the packet).  The coarse relaxation
over `rho in [0,1]` fails at `(4,2,6)`, where the value is `-65/56`, while
the actual value there is `45/3136`.  The smallest interior value on the
grid is `T(3,2;5) = 18`.  For large `N`, the positive term dominates by a
factor of about `N`.

## Status ledger (2026-09-27, main agent) -- SUPERSEDED by the updated ledger at the end

Every entry below was checked by the main agent.  Most were also checked by
the independent checker (FM-CHK).

**Infinite Q3 families proved.**  "Twist" means the sign-swapped family,
which also holds.

| family | words | reference |
|---|---|---|
| all labels in `{1,2}` | any | repo Cor. 24B17B / Prop. 24B3 |
| two minus labels, `<= 7` factors | any | repo Cor. 23A9ZZ10 |
| twisted-normal-form classes | classes 1 and 2 | FM11 |
| `{1,3}` with a single 3 | any | FM18 |
| minus `{q,1,1,1}` + plus 1s | any `q`, twist | FM19 |
| minus `{s+1,t+1,1,1}` + plus 1s | `min(s,t) <= 2` | FM26 |
| same-sign `{1,3}` | `#3 <= 6` or even | FM20 |
| same-sign `{1,3}` | 1s on the matching side dominate | FM22 |
| `{1,2,3}` with 3s of one sign | `#1^+ >= #1^- + #3^-`, twist | FM28 (contains FM22) |
| `{1,2,3}` with one 3 | `<= 3` two's (any 1s), or hard half with odd `#2 >= 2r+2` | FM25, FM29 |
| `{1,2,3}` with one plus 3 and one minus 3 | `<= 3` two's, and further ranges | FM25, FM29 |
| `[1^-,1^-,j^+,1^+ x a]` | all `j`, `a` | FM24 |
| one arbitrary label at any `r` | label `>= (N-2)/3`, or `a ~ t` | FM27 |
| `[1^-,1^-,p^+,q^+,1^+ x a]` | `p+q >= a+2`, `p = q`, `p` or `q = 1`, `a = 0` | FM30 |
| `[1^-,1^-,1^+,2^+ x b,3^+]` | any `b` | FM21 |

**Mechanism.**

- *Proved:* FM28, the recursion cone for labels `<= 3` on a half-quadrant.
- *Conjectured:* Conjecture HQ, which extends the FM28 cone to every
  generator except the plus-odd labels on `m >= r` (and the twist).  Exact
  tests are clean so far.
- *Killed:* K15 (orthant dominance), K16–K19 (certificate cones,
  Rodrigues rotation, orbit grouping, PF).

**Open targets and assigned agents.**

| target | agent |
|---|---|
| proof of HQ | venus |
| adversarial test of HQ, then its all-minus part | saturn |
| the hard half, outside both cones | neptune |
| one-label Turán middle range | mars |
| two-label strict interior | jupiter |
| checks | mercury |

**Uncovered by every current route:** superseded; see the updated ledger at
the end of the note.  The earlier claim was false as written: FM39 covers,
for example, `[3^+, 5^-, 3^-]`.  T1 and T2 in
general remain open.

### Proposition FM31 (the second Sp(4): the `u^2 v^2` tilt is Sp(4) Haar; main agent)

Put `alpha = theta+phi` and `beta = theta-phi`, so that `s = cos alpha` and
`t = cos beta`.  Then

```text
(x^2 - y^2)^2 = 16 sin^2 alpha sin^2 beta,
sin^2 theta sin^2 phi = (cos alpha - cos beta)^2/4.
```

So the doubled-SU(2) Haar measure tilted by `u^2 v^2 = S_1^2 D_1^2` is the
Sp(4) Weyl density in the angles `(alpha, beta)`.  Call this group `Sp(4)'`.
In `Sp(4)'`, with irreducibles `V'_(a,b)`, `a >= b >= 0`:

```text
2T = chi(C^4) = V'_10,      S = S_2 = chi(Lambda^2 C^4) = V'_11 + V'_00,
u^2 = det(1+g') = V'_11 + 2V'_10 + 3,   v^2 = det(1-g') = V'_11 - 2V'_10 + 3.
```

Hence, for `m, r >= 1`,

```text
Q3 value  ∝  < Ptilde det(1+g')^(m-1) det(1-g')^(r-1), 1 >_(Sp(4)'),
```

where `Ptilde` is the non-1 part of the word.  The generators decompose as
follows (exact decomposition by the Weyl formula, labels `<= 9`):

```text
H_(2k)/u    = V'_(k-1,k-1)                                   (genuine, irreducible)
H_(2k+1)    = V'_(k,k) + V'_(k,k-1) + V'_(k-1,k-1)            (genuine)
S_(2k)      = V'_(k,k) - V'_(k,k-2) + V'_(k-1,k-1) + V'_(k-2,k-2)   (one negative term; S_2 = V'_11 + V'_00)
S_(2k+1)/u  = V'_(k,k) - V'_(k,k-1) + V'_(k-1,k-1)                   (one negative term)
```

**Consequences.**

- *Every minus generator is a genuine `Sp(4)'` character.*
- *FM28 at `m = r = 1` is trivial*, because `T` and `S` are genuine
  characters.  For general `m >= r` the recursion handles the virtual weight
  `det(1-g')^(r-1)`.
- *Words with exactly two minus labels and only plus 1s* are trivially
  `>= 0`.  This agrees with FM11.
- *Conjecture HQ at `r = 1`* becomes a multiplicity-domination statement.
  For `G` in the cone generated by the `H`'s, `V'_10`, `V'_11` and
  `det(1+g')`, each plus-even label `S_(2k)` needs the pairing with
  `V'_(k,k) + V'_(k-1,k-1) + V'_(k-2,k-2)` to dominate the pairing with
  `V'_(k,k-2)`.  This is the same kind of statement as T1 (`m11 <= m20`),
  now in `Sp(4)'`.  The plus-odd labels carry `-V'_(k,k-1)`, a spin-type
  weight.  This is presumably why they fail on `m >= r`.

Via `Sp(4)' ~ Spin(5)`, the diagonal irreducibles `V'_(k,k)` are the SO(5)
spherical harmonics of degree `k`.

Open question: is there a 3-nesting-free (Pfaffian) or crystal model in
`Sp(4)'` in which the `V'_(k,k-2)` domination holds termwise?

### Corollary FM32 (more of {1,2,3} and a first label-5 case; FM-L123c, luna_max_neptune; accepted)

**Proved,** each by an explicit closed form with an all-positive numerator
certificate:

- `E[S^4 G] >= 0` for all `m, r`.  With FM29, one plus 3 is proved for up
  to four 2s, for all 1s.
- `E[S^4 G H_3] >= 0` for `m, r >= 1`.  The mixed family is proved for up
  to four 2s.
- `E[Q_5] >= 0` and `E[S Q_5] >= 0` for `m > r`, where `Q_5 = S_5/u`.  This
  is the first label-5 plus-odd case on the hard half.

**Checks.**  The main agent verified all four closed forms exactly on grids
up to 6 or 7.

**Next open formulas.**

- `E[S^5 G] >= 0` for `m > r >= 2`;
- `E[S^5 G H_3]` for `N >= 5`;
- `E[S^2 Q_5]` and `E[Q_7]` on `m > r`.

Growth by exponent: each step here is one more exponent.

### Lemma FM33 (Turán middle range, partial; FM-T1-AB3, luna_max_mars; accepted)

In the middle range `0 < m = 2k-N < (N-2)/3`, with `x = c_k` and
`y = c_(k+1)`:

```text
D_k - D_(k+1) = ((m+2)/(k+2)) (x - (m+1)d y/(2(m+2)(N-k+1)))^2
              + (m/(N-k+1)) (1 - (m+1)^2 d^2 / (m(m+2)(N-m+2)(N+m+4))) y^2,
d = a - t.
```

- **Case `min(a,t) <= 3`: proved.**  There are exact closed forms for
  `(D_k - D_(k+1))/c_k^2`, all positive.  The main agent checked them in
  1,001 exact cases.
- **Case `min(a,t) >= 4`: proved in an explicit band,**
  `v_- <= (m+1)^2 <= v_+` with `v_+- = 11N - 27 +- 2 sqrt(30(N-2)(N-3))`.

**Consequence.**  The one-arbitrary-label family, at every `r`, is now
proved whenever:

- the number of 1s is at most 3, or
- `t <= 3`, or
- the label is `>= (N-2)/3`, or
- the parameters lie in the band.

**Open:** `min(a,t) >= 4`, outside the band and below `(N-2)/3`.  This is
exactly mars's residual inequality.

### FM-MECH-A2 (luna_max_saturn, adversary): a misread HQ is false; HQ as intended stands

- Saturn's counterexample `E_(2,2)[T D_5^2] = -1/2` multiplies the `v^2` of
  the two minus 5s into the integrand at `(m,r) = (2,2)`.  In total-exponent
  terms this is `(m,r) = (2,3)`, with `m < r`, where `T`-odd terms carry the
  sign of `Delta`.  Likewise `E_(1,1)[T D_1^2] = -1` is the point `(1,2)`.
- Under the intended convention, the same product gives
  `E_(2,2)[T Phi_5^2] = +1` (main agent, exact).
- Saturn's screen of the one-atom layer (13,056 cases) is clean.  Its
  two-atom negatives (651) use the misread convention and must be re-run.
- Useful consequence: the cone is **not** closed under multiplication by the
  prefactors `u`, `v`.  Only the total-exponent formulation is a candidate.

### FM-CHK6 (luna_max_mercury, independent checker) and repairs

- **FM23 (Newton step): ACCEPTED.**  The checker gave its own proof via the
  derivative `P^(k-1)` and Cauchy--Schwarz.  The earlier FM-CHK5 dispute is
  resolved in the main agent's favour.  Newton gives `Delta_m >= 0`; the
  monotonicity is the separate FM27/FM33 question.
- **FM24: ACCEPTED.**
  - Explicit statement for odd labels: `T_1(j;a) = Delta_(j/2) - Delta_(j/2+1)`
    holds on the half-integer grid `j/2 in Z + 1/2` when `a` is odd.  The
    value `j = 0` uses `Shat_0 := 2`.
  - `T(p,q;0) = 3[p=q] + [|p-q|=2] - 2[p=q=1]`.
- **FM25: ACCEPTED.**  The closed forms, with `d = m-r >= 0` and
  `N = d+2r` (FM25 used `Delta` for `d`):

  ```text
  E[S G]   = 4 P_1(d,r)/((N+2)(N+3)^2(N+4)),
  P_1 = 3d^4 + 4d^3 r + 6d^3 + 16d^2 r^2 + 32d^2 r + 9d^2 + 16d r^3 + 24d r^2
        + 4d r + 6d + 16r^4 + 64r^3 + 100r^2 + 84r + 36;
  E[S^2 G] = 8 P_2(d,r)/((N+2)(N+3)^2(N+4)^2(N+5))   (P_2 as in FM25's packet, all coefficients positive);
  E[G H_3] = 4(4D^4 + 6D^2 N + 2D^2 + N^4 + 6N^3 + 11N^2 + 12N + 18)/((N+2)(N+3)^2(N+4)),
  E[S G H_3] = 8(8D^6 + 4D^4 N^2 + 72D^4 N + 136D^4 + 2D^2 N^4 + 26D^2 N^3 + 174D^2 N^2
             + 502D^2 N + 528D^2 + N^6 + 13N^5 + 71N^4 + 233N^3 + 480N^2 + 450N)
             /((N+2)(N+3)^2(N+4)^2(N+5)),   D = m-r.
  ```

  Word counts, with `a = #2^+`:

  | 3s in the word | `2m` | `2r` |
  |---|---|---|
  | one `3^+` | `#1^+ + #2^- + 1` | `#1^- + #2^-` |
  | one `3^-` | `#1^+ + #2^-` | `#1^- + #2^- + 1` |
  | one of each | `#1^+ + #2^- + 1` | `#1^- + #2^- + 1` |

- **FM26: `t = 1` ACCEPTED.**  It follows from FM19 by a binomial ratio.
  For `t = 2` the denominator display is now complete, and `Pi_2` is as
  displayed in the FM26 section.  A symbolic check of `t = 2` by the
  checker is pending; the main agent checked it exactly on the grid
  `s <= 11`, `v <= 7`.
- **M0: scalar `2^n (-1)^d` restored.**  The PSD claim is now stated for
  even `d` only.
- **M1: ACCEPTED.**  Boundary positivity is at `x = 2` or `y = 2`; at
  `x = -2` the generators carry the sign `(-1)^p`.
- **K15: confirmed with an exact witness.**
  `int_(u,v >= 0) S_3 dsigma dsigma = -256/(1575 pi^2) < 0`.

### Corollary FM34 (absorption at `r = 1`; FM-L123d (neptune) + main agent)

**Neptune (proved).**  In `Sp(4)'`, with `P_3 = S_3/u`:

```text
P_3 H_3 = (V'_11 + V'_00)^2 - (V'_10)^2 = V'_22 + V'_11 + V'_00,
```

which is genuine.  Hence, for `r = 1` and all `m >= 1`, `a >= 0`,

```text
P_3 H_3 S^a det(1+g')^(m-1)
```

is genuine.  So Q3 holds for every word with one plus 3, one minus 3, any
number of plus 2s and plus 1s, and exactly one further unit of `v`
(`#1^- + #2^- = 1`).

**Neptune (killed).**  Multiplicity domination for an *arbitrary* genuine
`G` is false: `<A_1 V'_20 det(1+g'), 1> = -1`.  Conjecture HQ restricts `G`
to the word cone, so it is not affected.

**Absorption table (main agent, exact in `Sp(4)'`, labels `<= 9`).**

- `S_p H_q` (with the `u`-quotients) is a **genuine** `Sp(4)'` character iff
  `q >= p`, and in addition `q` is odd when `p` is odd.
- `S_2` is genuine by itself.  No product of two plus labels `>= 3` is
  genuine.
- **Reason, in the original Sp(4).**  By FM5,
  `S_p D_q = sum_(j=|p-q|)^(p+q) D_j + (V_q x V_p - V_p x V_q)`.
  The bracket equals `D_1 Res V_(q-1,p)` when `q > p`, and 0 when `q = p`.
  So `S_p H_q = sum_j H_j + Res V_(q-1,p)`, which is genuine in the original
  Sp(4) for `q >= p`.

**Consequence (`r = 1`, two minus labels `q_1, q_2`).**  Q3 holds whenever
every plus label `>= 3` can be matched injectively to a minus label
`q >= p` of the right parity (at most two such labels), and in addition:

- *In the original Sp(4):* all other plus labels are 1s.
- *In `Sp(4)'`:* all other plus labels are 1s or 2s, `m >= 1`, and plus
  4s, 6s, ... must also be matched.

The general statement needs a proof for all labels (assigned).

**Obstruction for `r >= 2`.**  The extra weight `v^2 = det(1-g')` is never
absorbed.  `det(1-g') H_q H_(q')` is non-genuine for every `2 <= q, q' <= 9`.
`det(1-g') det(1+g') = det(1-g'^2)` is non-genuine too, so
`E_(2,2)[G] >= 0` fails for some genuine `G`.  Beyond `r = 1`, positivity
needs the recursion mechanism (FM28, HQ), not genuineness.

### Kill K20 (Conjecture HQ as a cone with auxiliary `T` is false; FM-MECH-A3, luna_max_saturn)

The counterexample is at total exponents `(m,r) = (3,3)` (validity holds,
`m = r`):

```text
E_(3,3)[T^3 Phi_11] = -1/400,
Phi_11 = S^5 + 2S^4 T - 16S^3 T^2 - 4S^3 T - 4S^3 - 24S^2 T^3 + 24S^2 T^2
       + 48S T^4 + 64S T^3 + 3S + 32T^5 - 64T^4 - 64T^3 + 6T.
```

Here `Phi_11` is the minus-11 quotient.

- The main agent confirmed the value exactly.
- The one-atom screen is clean for every label `<= 10` (28,160 cases).  The
  earlier screens stopped at label 8, so they could not see this.
- `T^1 Phi_11` and `T^2 Phi_11` stay `>= 0` for `r <= m <= 15`.
- Word-level expressions without auxiliary `T`, namely `E[S^b Phi_q]` for
  `q <= 15`, `b <= 3`, `r <= m <= 12`, stay `>= 0` (main agent).

**Status.**

- The recursion-cone mechanism is valid for labels `<= 3` (FM28) and does
  **not** extend with the monomial cone `T^a S^b`.
- The word-level statement: Q3 for words with no plus-odd label `>= 3` and
  `#1^+ >= #1^- + #(odd minus labels)`.  It is just Q3 on that class, with
  no mechanism attached, and it is not contradicted.
- The mechanism picture is now:
  - the recursion cone for small labels;
  - genuineness with absorption at `r = 1` (FM31, FM34);
  - no mechanism yet for `r >= 2` with large labels.

### FM-T1-AB4 (luna_max_mars): one-label residual region pinned down

After FM27 and FM33, the only open case of the one-arbitrary-label family
(every `r`) is

```text
a, t >= 4,   0 < m = 2k-N < (N-2)/3,   (m+1)^2 (a-t)^2 > m(m+2)(N-m+2)(N+m+4),   c_(k+1) != 0,
```

where the target is the ratio inequality stated in the packet.  The region
is thin:

- for `m = 1` it requires roughly `|a-t| > 0.866 N`, so `t < 0.067 N`, and
  first occurs at `N = 79` (`a = 75`, `t = 4`, exact margin `55/7776`);
- for `m = 2` it requires `t < 0.028 N`, first at `N = 190`.

An exact scan over `4 <= a <= 500`, `4 <= t <= 21` (96,929 residual cases)
found every case positive.  The minimum of `Delta_k/(|x|+|y|)^2` is
`7/127755`, at `(a,t,k) = (500,5,253)`.  Near the centre this is the
Krawtchouk-to-Hermite regime (`t << N`).

### FM-MECH-B3 (luna_max_venus): the r=1 plus-even domination, no proof

- **The `r = 1` reduction is confirmed:**
  `<G S_(2k), 1> = m_G(V'_kk) + m_G(V'_(k-1,k-1)) + m_G(V'_(k-2,k-2)) - m_G(V'_(k,k-2))`.
- **Exact screens are clean.**
  - A single `S_(2k)`, `k = 2..6`, against `G` a product of at most four
    generators: 3,575 margins, all `>= 0`.  The minimum margins are
    `1,1,1,2,2`.
  - Two plus-even labels: 1,638 cases, minimum 0.
- **Two shortcuts killed.**
  - `S_4 V'_10 = V'_32 - V'_30 + V'_21 + V'_10` is not genuine.
  - `S_4^2 = V'_44 - V'_42 + 2V'_40 + 2V'_33 - V'_31 + 3V'_22 - 3V'_20 + 3V'_11 + 4V'_00`
    is not genuine, although `<S_4^2, 1> = 4`.

**Core for `r >= 2` (main agent).**  In `Sp(4)'`, each further pair of
minus 1s multiplies by `det(1-g') = 3 + V'_11 - 2V'_10`.  For words with
exactly four minus labels and plus labels 1, Q3 becomes
`m'_11 + 3m'_00 >= 2m'_10` for `M = prod_(j<=4) H~_(q_j) det(1+g')^(m-1)`.
In the original Sp(4) it is `5m_00 - 3m_11 + m_20 >= 0` for
`Sym^(k_1) (x) ... (x) Sym^(k_4) (x) (C^4)^(x)a`.  Assigned as FM-T1-4S.

### FM30 extension (FM-T2-PQ4, luna_max_jupiter; accepted)

In the strict interior, put `g = p-q` and `h = (a+2-p-q)/2 >= 1`.  Then
`sgn M = sgn Pi`, where

```text
Pi = 12h^2 - 4A h + C,
A  = g^2 + 2gq - g + 2q^2 - 2q - 3,
C  = (g+2q+1)(g+2q+2)(g^2-g-2q).
```

**Proved: every interior case with `Pi >= 0`.**  Here `M >= 0`, and
`Delta_d - Delta_(e+1) = sum_j T_1(j;a) >= 0` by FM24.  The packet checked
the normalisation against the ballot formula on 378 exact triples; the
smallest value is `T(3,2;5) = 18`.

**Open.**  For each `(g,q)` there is a bounded window of `h`, the roots of
`Pi`, where the one-variable inequality
`rho <= (sqrt(gamma^2 + 4 delta_d delta_(e+1)) + gamma)/(2 delta_(e+1))`
remains.  For example, `(p,q) = (3,2)` leaves `h = 1..4`.

### Corollary FM35 (two symmetric powers, smaller one equal to 3; FM-T1-2L3, luna_max_saturn; accepted)

Put `n = s-3` and `C_3 = binom(s-3+2v, s+v-3)`.  Then

```text
L(s,3,s-3+2v) = C_3^2 (s+1)(s+2)(s+3)(s+2v)(s+2v-1)(s+2v-2)^2 Pi_3(n,v)
   / ((s+v)^2(v+1)^2(v+2)(s+v-2)^2(s+v-1)^2(s+v+1)^2(s+v+2)^2(s+v+3)^2(s+v+4)),
Pi_3 = (n+2)(n+3)(n+4)(n+5)(n+6)(n+7)(n^2-3n+6)
     + (4n^7+63n^6+297n^5+395n^4+2083n^3+19318n^2+51760n+44496) v
     + (9n^6+105n^5+460n^4+2273n^3+11591n^2+24898n+20520) v^2
     + 2(5n^5+20n^4+3n^3+886n^2+4198n+5286) v^3
     + v^4 (A(n) + B(n) v + 300 v^2),
A = 5n^4-50n^3+103n^2+3230n+5388,   B = -60(n^2-5n-23),
B^2 - 1200A = -1200(2n^2((n-5)^2+58) + 2540n + 3801) < 0.
```

The main agent checked the closed form exactly for `3 <= s <= 12`,
`v <= 6` (70 cases).  The real roots of the discriminant polynomial are
negative.

**Corollary FM35.**  Q3 holds for minus labels `{s+1, 4, 1, 1}` with any
number of plus 1s.  With FM19 and FM26, the case `min(s,t) <= 3` is
complete.

The slices `t = 1, 2, 3` share visible structure: `C_t^2` times linear
factors in `s+2v`, times a polynomial `Pi_t` that is positive by a
discriminant argument.  A general-`t` form is assigned as FM-T1-2L4.

### FM-T1-4S (luna_max_venus): the pointwise-square family `{q,q,2,1}`

Minus labels `{q,q,2,1}` with any `a` plus 1s give

```text
J = int D_q^2 (x-y)^2 (x+y)^(a+1) dsigma dsigma.
```

This is `>= 0` pointwise for odd `a`, and 0 for even `a`.  It is correct but
elementary, in the same class as the equal-pairs case of FM-T1-SX.  Example:
`J(4,1) = 8`.  The general four-minus core is untouched.

### Lemma FM36 (one-label family: labels 1, 2, 3 in full; FM-T1-AB5, luna_max_mars; accepted)

At the target index `k = (N+m)/2`, reciprocity (`c_(N-j) = (-1)^t c_j`)
together with the recurrence determines every coefficient near the centre
from a single value.  Hence `D_k - D_(k+1) = Q_m(N,t) c^2`, with (notation repaired after FM-CHK8):

```text
m=1:  Q = 32(t+1)(N-t+2)/((N+3)^2(N+5))            [t even],   32(t+2)(N-t+1)/((N+3)^2(N+5))  [t odd],  c = c_k;
m=2:  Q = 16(t+1)(N-t+1)(N^2+6N+8(t-N/2)^2)/((N+2)^2(N+4)^2(N+6))   [t even, c = c_(N/2)],
      Q = 32(t+2)(N-t+2)/((N+4)^2(N+6))                              [t odd,  c = c_(N/2+1)];
m=3:  Q = 64(t+1)(N-t+2)((N-2t)^2+4t^2+12t-1)/((N+3)^2(N+5)^2(N+7))            [t even],
      Q = 64(t+2)(N-t+1)(5N^2-12Nt+12N+8t^2-12t-1)/((N+3)^2(N+5)^2(N+7))       [t odd],  c = c_((N+1)/2).
```

All of these are `>= 0`; for `m = 3` this uses `t >= 4`, and FM33 covers
`t <= 3`.  The main agent checked all three exactly: 3,729 cases,
`4 <= N < 70`, all `t`.

**Consequence.**  In the one-arbitrary-label family, at every `r`, labels
`q <= 3` (and `s+1 <= 3` for Sym) are proved for all `a` and `t`.

**Open.**  Labels `q = m >= 4` with `min(a,t) >= 4`, `q < (N-2)/3` and
`F_q(N,l) = (q+1)^2(N-2l)^2 - q(q+2)(N-q+2)(N+q+4) > 0`.  Correction after
FM-T1-AB8: this does **not** force `N` large.  `F > 0` holds iff
`N < N_-` or `N > N_+`, where
`N_+- = 2vl + 3v - 3 +- 2 sqrt(v(v-1)(l+1)(l+2))` and `v = (q+1)^2`.  The lower
branch is non-empty: `(m,N,t) = (65,199,4)` lies in it, with exact margin
`279071056950/4990837309097 > 0`.  The first instance is `(a,t,N,k,m) = (538,4,542,273,4)`,
exact margin `1816764288/4713370019435 > 0`.

**Main-agent observation on general `m` (for FM-T1-AB6).**  Work symbolically,
with `N` even, `t` even and `d = N-2t`.  Then

```text
Delta_((N+m)/2) - Delta_((N+m)/2+1) = 16(t+1)(N-t+1) P_m(N,d) c_(N/2)^2
                / ((N+2)^2 (N+4)^2 ... (N+m+2)^2 (N+m+4)),
P_2 = N^2 + 6N + 2d^2,
P_4 = 3N^4 + 36N^3 - 6N^2 d^2 + 108N^2 - 36N d^2 + 8d^4 - 128d^2,
P_6 = 3N^6 + 54N^5 + 12N^4 d^2 + ... + 32d^6,   P_8 = 5N^8 + ... + 128d^8.
```

The top homogeneous part is `N^m (sum_(j=0)^(m/2) T_(2j)(x) + m/2)` with
`x = d/N`, where `T` are the Chebyshev polynomials.  This was checked for
`m = 2, 4, 6, 8`.

- Since `sum_(j=0)^M cos(2j theta) = (1 + sin((2M+1)theta)/sin theta)/2 >= -0.43 M`,
  the top part is `>= ~0.57 (m/2) N^m > 0`.
- A proof for all `m` then needs control of the lower-order terms in the
  residual regime `N >~ 4(m+1)^2 l`.  Alternatively, the full `P_m` may
  have an exact Chebyshev/Dirichlet form.

### Lemma FM37 (one-label family: labels 4 and 5; FM-T1-AB6, luna_max_mars; accepted)

**General central recursion.**  `D_k - D_(k+1) = Q_m b^2` (repaired after FM-CHK8), where

```text
Q_m = R_r^2 - R_(r-1) R_(r+1) - R_(r+1)^2 + R_r R_(r+2),
R_(j+1) = (d R_j - (n+eta-j+1) R_(j-1))/(n+j+1).
```

Here `N = 2n + eta` and `r = (m+eta)/2`.  The base and initial ratios:

- `N` even, `t` even: `b = c_n`, `R_0 = 1`, `R_1 = (N-2t)/(N+2)`;
- `N` even, `t` odd: `b = c_(n+1)`, `R_0 = 0`, `R_1 = 1`;
- `N` odd: `b = c_n`, `R_0 = 1`, `R_1 = (-1)^t`.

**Closed forms for `m = 4, 5`** (`t = min(a,t) >= 4`, `d = N-2t`):

- `m = 4`:
  - `t` even: `P_4e = 5d^4 + d^2(48t^2+144t-20) + 96dt^3 + 432dt^2 + 432dt + 48t^4 + 288t^3 + 432t^2`;
  - `t` odd: `P_4o = 5d^2 + 4dt + 6d + 4t^2 + 12t - 16`.
- `m = 5`:
  - `t` even: `P_5e = 3d^4 - (8t+12)d^3 + (72t^2+216t+6)d^2 + (96t^3+432t^2+440t+12)d + 48t^4+288t^3+216t^2-648t-9`.
    Its `d`-quadratic `3d^2-(8t+12)d+(72t^2+216t+6)` has discriminant `-800t^2-2400t+72 < 0`.
  - `t` odd: `P_5o = 35d^4 + (56t+84)d^3 + (72t^2+216t-218)d^2 + (96t^3+432t^2-8t-660)d + 48t^4+288t^3+216t^2-648t-9`.
    All of its coefficients are positive for `t >= 4`.

With prefactors, `Q = (D_k - D_(k+1))/b^2`:

```text
Q_4e = 16(t+1)(d+t+1) P_4e / ((N+2)^2(N+4)^2(N+6)^2(N+8)),
Q_4o = 32(t+2)(d+t+2) P_4o / ((N+4)^2(N+6)^2(N+8)),
Q_5e = 32(t+1)(d+t+2) P_5e / ((N+3)^2(N+5)^2(N+7)^2(N+9)),
Q_5o = 32(t+2)(d+t+1) P_5o / ((N+3)^2(N+5)^2(N+7)^2(N+9)).
```

All four are positive for `t >= 4`.  The main agent checked them exactly:
1,722 cases, `N < 90`.

**Consequence.**  In the one-arbitrary-label family at every `r`, labels
`<= 5` are proved for all `a`, `t`.

**Open.**  `m >= 6`, `t >= 4`, `m < (N-2)/3`, `F_m > 0`.  The first
instance is `N = 1068`, `(a,t,m) = (1064,4,6)`, with a positive margin.

### FM-CHK7 (luna_max_mercury) and repairs

- **FM31: ACCEPTED for all `k`.**  The checker proved the decompositions
  with generating functions:

  ```text
  sum_k H_(2k+1) q^k = Q/D,
  sum_k (H_(2k+2)/u) q^k = (1+q)/D,
  Q = 1 + (2+2T)q + q^2,
  D = Q^2 - u^2 q (1+q)^2.
  ```

  These match the Weyl-alternant sums; the `S_p` forms then follow via
  FM6.  The normalisation is `dmu_(Sp(4)') = (1/2) u^2 v^2 dsigma dsigma`,
  so `int Ptilde u^(2m) v^(2r) = 2 <Ptilde det(1+g')^(m-1) det(1-g')^(r-1), 1>`.
- **FM28: ACCEPTED.**  The cone-membership table is confirmed through
  label 7.
- **FM27: ACCEPTED.**  This covers (FM27a), (FM27b) via the weights
  `lambda_j = 1/binom(N+1,j+1)`, and both reductions.
- **FM33: square identity and band ACCEPTED.**  Repair for
  `min(a,t) <= 3`: `c_k` can vanish (e.g. `a = 2`, `t = 14`, `k = 10`).  In
  that case the square identity gives `Delta_k = m c_(k+1)^2/(N-k+1) >= 0`
  directly.  When `c_k != 0`, the ratio forms from the FM-T1-AB3 packet
  apply:

  ```text
  l=0: Delta/c_k^2 = 16(N+1)(N+2)(m+1)/((N-m+2)(N+m+2)^2(N+m+4))
  l=1: 16N(N+1)(m+1)(m+2)/(m(N-m+2)(N+m+2)^2(N+m+4))
  l=2: 16N(N-1)(m+1)A_2/((N-m+2)(N-m^2)^2(N+m+2)^2(N+m+4)),
       A_2 = (m^2+2m-N)^2 + 2N^2 + 6N - 4m^2 - 8m
  l=3: 16(N-2)(N-1)(m+1)(m+2)A_3/(m(3N-m^2-2)^2(N-m+2)(N+m+2)^2(N+m+4)),
       A_3 = (m^2-3N)^2 + 6N^2 - 12Nm + 4m^3 - 4m^2 - 16m + 18N
  ```

  These were checked in 1,001 exact cases.
- **FM30: the closed form for `M` is now a proved symbolic identity**
  (main agent, sympy).  The 12-term reflection stencil, divided by
  `beta_(e-1) beta_(d-1)` through exact binomial shift ratios, equals the
  displayed rational function identically in `(N,d,e)`.  With the boundary
  and `Pi >= 0` proofs, this repair is complete.
- **FM26 (`t = 1, 2`) and FM35 (`t = 3`): proved symbolic identities**
  (main agent, sympy).  The sum over constituents of the 24-term `Lambda`
  formula, divided by `C_t^2` through exact binomial shift ratios, equals
  the displayed closed forms identically in `(s,v)`.  Out-of-range
  binomials vanish correctly inside the ratio products, and the
  denominators are positive on the domain.

### Corollary FM38 (two symmetric powers, smaller one at most 5; FM-T1-2L4, luna_max_saturn; accepted)

Put `n = s-t`, `a = n+2v` and `C_t = binom(n+2v, n+v)`.  For `t = 4, 5`:

```text
L(s,t,a) = C_t^2 (s+1)(s+2)(s+3) F_t(s,v) Pi_t(n,v) / D_t(s,v),
F_4 = (s+2v-3)(s+2v-2),   F_5 = (s+2v-3)(s+2v-2)(s+2v-4)^2,
D_t = (s+v)^2 (v+1)^2 (v+2) prod_(h=1)^(t-1) (s+v-h)^2 prod_(h=1)^3 (s+v+h)^2 (s+v+4),
Pi_t = sum_(j=0)^10 A_(t,j)(n) v^j.
```

The coefficients (inserted after FM-CHK8) are as follows.  For `t = 4`
(`n = s-4`):

```text
A_(4,0)  = (n+1)(n+2)(n+3)^2(n+4)^2(n+5)(n+6)(n+7)(n+8)(n^2-3n+6)
A_(4,1)  = 2(n+3)(n+4)(5n^9+137n^8+1467n^7+7766n^6+22545n^5+54689n^4+204595n^3+606216n^2+833172n+403056)
A_(4,2)  = (n+3)(46n^9+1194n^8+12217n^7+64055n^6+211125n^5+712255n^4+2856700n^3+7934928n^2+11086344n+5834016)
A_(4,3)  = 2(64n^9+1510n^8+14154n^7+71674n^6+272911n^5+1185936n^4+4877611n^3+12859144n^2+17850276n+9854520)
A_(4,4)  = 239n^8+4402n^7+31942n^6+147414n^5+737342n^4+3610862n^3+11445217n^2+19217406n+13013016
A_(4,5)  = 2(154n^7+2016n^6+11393n^5+69476n^4+438237n^3+1648378n^2+3198378n+2542878)
A_(4,6)  = 2(133n^6+868n^5+4500n^4+62066n^3+344267n^2+779190n+661266)
A_(4,7)  = 28(5n^5-20n^4+121n^3+3542n^2+12012n+12690)
A_(4,8)  = 7(n+3)(5n^3-165n^2+1366n+4464)
A_(4,9)  = -420(n^2-13n-31),   A_(4,10) = 2100
```

For `t = 5` (`n = s-5`):

```text
A_(5,0)  = (n+2)(n+3)(n+4)^2(n+5)^2(n+6)(n+7)(n+8)(n+9)(n^2-3n+6)
A_(5,1)  = (n+4)(n+5)(8n^9+263n^8+3402n^7+21712n^6+72202n^5+167829n^4+709860n^3+2917236n^2+5550432n+3879936)
A_(5,2)  = 31n^10+1083n^9+15629n^8+122126n^7+597261n^6+2324683n^5+9553739n^4+35740732n^3+87617884n^2+115215632n+62659200
A_(5,3)  = 2(37n^9+1040n^8+11616n^7+68684n^6+286039n^5+1349418n^4+6571120n^3+20930066n^2+34801020n+22871520)
A_(5,4)  = 121n^8+2632n^7+22274n^6+113004n^5+616387n^4+3671048n^3+14447090n^2+29250084n+22912080
A_(5,5)  = 4(35n^7+546n^6+3553n^5+21995n^4+153486n^3+709213n^2+1745832n+1786590)
A_(5,6)  = 4(28n^6+245n^5+1560n^4+18739n^3+115016n^2+327012n+365490)
A_(5,7)  = 56(n^5-2n^4+35n^3+776n^2+2610n+2610)
A_(5,8)  = 14(n^4-26n^3+191n^2+1874n+4080)
A_(5,9)  = -168(n^2-11n-15),   A_(5,10) = 840
```

**Sturm data (main agent, sympy `count_roots(0, oo)`).**

- Every `A_(t,j)` with `j <= 8` has 0 roots in `[0, inf)`.  Its value at
  `n = 0`:
  - `t = 4`: `2903040, 9673344, 17502048, 19709040, 13013016, 5085756, 1322532, 355320, 93744`;
  - `t = 5`: `43545600, 77598720, 62659200, 45743040, 22912080, 7146360, 1461960, 146160, 57120`.
- `A_(t,9)^2 - 4 A_(t,8) A_(t,10)` has 0 roots in `[0, inf)`.  Its value
  at 0 is `-617929200` (`t = 4`) and `-185572800` (`t = 5`).
- The checker's independent values `L(4,4,0) = 3`, `L(5,4,1) = 2`,
  `L(5,4,3) = 6` are consistent with this.

**Checks (main agent).**

- Both closed forms are exact symbolic identities in `(s,v)`: the constituent
  sums of the 24-term formula, divided by `C_t^2` through binomial shift
  ratios.
- Sturm counts show that `A_(t,j)` has no root in `[0, inf)` for
  `j <= 8`, and that the discriminant of the final quadratic
  `A_(t,8) + A_(t,9) v + A_(t,10) v^2` has none either (its value at 0 is
  negative).  Hence `Pi_4 > 0` and `Pi_5 > 0`.

**Corollary FM38.**  Q3 holds for minus labels `{s+1, t+1, 1, 1}` with any
number of plus 1s whenever `min(s,t) <= 5`.

**Pattern for `t = 1..5`:**

- The denominator is `D_t`.
- `deg_v Pi_t` is `2, 6, 6, 10, 10`, with leading coefficients
  `60, 600, 300, 2100, 840`.
- No uniform-`t` certificate is known yet.  The case `t >= 6` is open.

### FM-T2-PQ5 (luna_max_jupiter): the two-label window, localised

The negative window `Pi < 0` lies in `h < 2A < 2(p^2+q^2)`, and requires
`p < 3q`.  It is infinite: `p = q+1`, `h = 1` gives
`Pi = -8q^3 - 28q^2 - 12q + 24 < 0` for every `q`.

Kills:

- the bound `rho <= ((h+q)/(h+p+q+1))^(q+1)`, and its Jensen refinement,
  are not shown to suffice;
- pairing `M` with only the two endpoint one-label terms fails, e.g. at
  `(p,q,h) = (6,5,7)`.

Still open.

### Theorem FM39 (absorption at `r = 1`, all labels; FM-L123e, luna_max_neptune; accepted)

**Branching identities from Sp(4) to `Sp(4)'`.**  Put `c = (A+B)/2` and
`d = (A-B)/2` (or the analogues below).

- For `A + B` even:

  ```text
  Res V_(A,B) = sum_(i=0)^1 sum_(j=max(0,d-1))^(min(d, c-1+i)) V'_(c-1+i, j)   (truncation added after FM-CHK8).
  ```

- For `A + B` odd, with `c = (A+B-1)/2` and `d = (A-B-1)/2`:

  ```text
  Res V_(A,B) = u V'_(c,d).
  ```

**Absorption.**

- *Original Sp(4):*
  `S_p H_q = sum_(j=|p-q|, step 2)^(p+q) H_j +/- Res V_(...)`, which is
  genuine iff `q >= p`.
- *`Sp(4)'`, normalised factors:* `P_p M_q` is genuine iff either `p <= 2`,
  or `q >= p` with `p` even or `q` odd.
- *Converse:* the negative constituent survives in every other case, by a
  central-parity or off-diagonal argument.

**Theorem FM39 (`r = 1`, exactly two minus labels).**  Q3 holds if either
of the following matchings exists:

1. In `Sp(4)'` (needs `m >= 1`): every plus label `p >= 3` is matched
   injectively to a minus label `q >= p`, with `q` odd when `p` is odd.
   Plus 1s and plus 2s are arbitrary.
2. In the original Sp(4): every plus label `p >= 2` is matched injectively
   to a minus label `q >= p`.  The other plus labels are 1s.

This adds coverage of unbounded length beyond Cor. 23A9ZZ10 (`<= 7`
factors) and FM11.

**Kill (auxiliary `T` at `r >= 2`).**  `G = P_3 H_5`
(`= V'_00 + V'_11 + V'_22 + V'_30 + V'_33`, genuine) gives
`E_(2,2)[G T] = -1/4`.  Marked-moment positivity with an auxiliary `T` is
false; the word-level statement is not affected.

**Word-level core at `r = 2`, `m = 2`.**  Here `CD = det(1+g') det(1-g') =`
`V'_22 - 3V'_20 + 2V'_11 + 6V'_00`, so Q3 becomes

```text
3 m_G(V'_20) <= m_G(V'_22) + 2 m_G(V'_11) + 6 m_G(V'_00)
```

for word modules `G`.  This is a T1-type inequality in `Sp(4)'`.

### FM-T1-AB7 (luna_max_mars): the one-label family beyond label 5

No new label was closed.  Main-agent data for the next attempt (`N`, `t`
even, `d = N-2t`): the `d^0` coefficient of `P_m` factors completely:

```text
N(N+6);   3N^2(N+6)^2;   3N^2(N-4)(N+6)^2(N+10);
5N^2(N-4)^2(N+6)^2(N+10)^2;   5N^2(N-8)(N-4)^2(N+6)^2(N+10)^2(N+14),
```

for `m = 2, 4, ..., 10`.  In the basis `prod_(i<j)(d^2 - (2i)^2)` the
remaining coefficients do not factor cleanly.

### FM-T1-AB8 (luna_max_mars): the residual has two branches

The residual has two branches:

- *Upper branch*, `N > N_+ ~ 4 v l`: `N` is large relative to `m^2`, the
  regime of the Chebyshev asymptotics.
- *Lower branch*, `3m+2 < N < N_- ~ v/(4(l+1))`: here `m` is large relative
  to `l` and `k = (N+m)/2` lies just below the FM27 threshold `(2N-1)/3`.
  The negative middle term `(3k-2N+1)` in (FM27b) is then small.

A refined tail estimate in (FM27b) is the natural tool for this branch.
Exact checks at the first points of each branch are positive.

### Lemma FM40 (two symmetric powers at `r = 2` telescope; FM-T1-2L5, luna_max_saturn; accepted)

For `s >= t`:

```text
H_(s+1) H_(t+1) = sum_(j=0)^t H_(s+t-2j+1) + sum_(j=1)^t sum_(h=1)^j Res V_(s+t-2j+h, h).
```

Hence, when `s+t+a` is **even**, with `c_n = [z^n](1+z)^a(1-z)^3` and
`Delta_k = c_k^2 - c_(k-1) c_(k+1)` (the FM27 one-label sequence at `r = 2`),
the formula below holds.  When `s+t+a` is odd, `L = 0` by central parity.
(Parity clause added after FM-CHK9.)

```text
L(s,t,a) = Delta_((a+s-t+4)/2) - Delta_((a+s+t+6)/2) + M_2(s,t,a),
M_2 = sum_(j=1)^t sum_(h=1)^j Lambda(s+t-2j+h, h; a),
```

and `M_2` has an explicit binomial formula (in the packet).  This is the
`r = 2` analogue of jupiter's `r = 1` split `T = Delta_d - Delta_(e+1) + M`.
The main agent checked the identity exactly in 324 cases (`a < 12`,
`t <= 5`).

`M_2` can be negative: `M_2(2,1,1) = -2` while `L = 2`.

**Open (`t >= 6`, interior):**
`M_2 >= Delta_((a+s+t+6)/2) - Delta_((a+s-t+4)/2)`.

### FM-T2-PQ6 (luna_max_jupiter): two plus labels, `q = 2, 3` complete

- **`q = 2, 3` proved.**  The negative-`Pi` windows are finite lists (`p < 3q`),
  and all 37 cases were evaluated exactly from the ballot formula; all are
  positive.  So the two-plus-label family `[1^-,1^-,p^+,q^+,1^+ x a]` is
  proved whenever `min(p,q) <= 3`.
- **`q >= 4` (conditional).**  Jensen gives
  `rho <= B = ((2h+q)/(2h+2p+q+2))^(q+1)`.  A sufficient condition is
  `R B + U B^2 <= 1`, with `R = N(p+1)(q+1)(-Pi)/(h S D)` and
  `U = E(h+p+1)(h+q+1)/(h S D)`.  An exact screen for `q = 4..12` found a
  maximum of `0.386` (at `(q,p,h) = (4,5,5)`).  A uniform proof is open.

*Main-agent scan (floating point, for guidance only).*  Over all negative-`Pi`
windows with `4 <= q <= 120`, the maximum of `R B + U B^2` stays at
`0.386`, attained at `(q,p,h) = (4,5,5)`.  Larger `q` only lowers it.  This
supports the plan: an analytic bound for `q >= Q_0`, plus an exact finite
check below.


## Status ledger, updated (2026-09-27, after FM-CHK8) -- SUPERSEDED by the ledger after FM54 at the end

Every entry was checked by the main agent, and most also by an independent
checker.  "Twist" means the sign-swapped family, which also holds.

**Proved infinite Q3 families.**

- **Labels in `{1,2}`:** all words.  Repo Cor. 24B17B / Prop. 24B3.
- **Two minus labels (`r = 1`):**
  - `<= 7` factors: repo Cor. 23A9ZZ10.
  - All plus labels 1: FM11.
  - Absorption: every plus label `p >= 3` matched injectively to a minus
    label `q >= p`, with `q` odd when `p` is odd; plus 1s and 2s arbitrary
    (`Sp(4)'` form, `m >= 1`).  Or, in the original form, every plus label
    `p >= 2` matched and the rest are 1s.  FM39.
  - `[1^-,1^-,j^+,1^+ x a]` for all `j`, `a`: FM24.
  - `[1^-,1^-,p^+,q^+,1^+ x a]`: **all `p, q, a`** (Theorem FM42).
  - `[q^-,1^-,p^+,1^+ x a]`: **all `p, q, a`** (Theorem FM49).
  - `[1^-,1^-,p^+,q^+,s^+,1^+ x a]` (`p >= q >= s`): `p+q-s > a+2` (FM48);
    exact screen of the residual for labels `<= 20`, `a <= 30`.
  - `[1^-,1^-,1^+,2^+ x b,3^+]`: FM21.
- **Labels `{1,3}` and `{1,2,3}`:**
  - A single 3: FM18.
  - Same-sign 3s: **all counts, all 1s** (FM18, FM20, FM22, **FM47**).
    With any number of 2s, on the matching half-quadrant (FM28).
  - One 3 with at most four 2s, and one plus 3 with one minus 3 with at most
    four 2s: FM25, FM29, FM32.  Also the plus-3/minus-3 family at `r = 1`
    (FM34), and further ranges.
- **One arbitrary label at every `r`** (`Shat_q` or `Sym^s` with fundamental
  labels):
  - label `<= 5`, all `a`, `t`: FM36, FM37;
  - label `>= (N-2)/3`: FM27;
  - `min(a,t) <= 3`, or inside the FM33 band.
  - Open: label `>= 6` on the residual's two branches.
- **Four minus labels `{s+1, t+1, 1, 1}` with plus 1s:** `min(s,t) <= 11`
  (FM19, FM26, FM35, FM38, FM41, FM43, FM46), plus the square case `s = t`,
  `a = 0`, the support boundary and odd parity.  Open: `s > t >= 12`,
  `a >= 1`, `s-t < a+2`, `s+t+a` even.  The telescoping
  identity FM40 holds for all `t` in even parity; in odd parity `L = 0`.
- **Minus labels `{q,1,1,1}` with plus 1s:** FM19, and the twist.
- **Elementary pointwise families:** equal pairs, and `{q,q,2,1}`.

**Mechanisms.**

- *Proved:*
  - the recursion cone for labels `<= 3` on a half-quadrant (FM28);
  - genuineness plus absorption at `r = 1` in `Sp(4)'` (FM31, FM34, FM39);
  - the Turán/telescoping structure for one and two labels (FM23, FM27,
    FM30, FM40).
- *Killed:* K15 (orthant dominance), K16–K19 (certificate cones), K20 (HQ
  with auxiliary `T`), marked-`T` positivity at `r >= 2` (FM39).

**Open core.**

- `r >= 2` with large labels.  In `Sp(4)'` this is
  `<Ptilde det(1+g')^(m-1) det(1-g')^(r-1), 1> >= 0`.  At `m = r = 2` it
  becomes `3 m20' <= m22' + 2 m11' + 6 m00'`, proved so far only for
  products of two even minus generators (FM-L123f).
- T1 and T2 in general.

### Round notes (after FM-CHK8)

- **FM-T1-AB9 (mars).**  The two-tail quadratic from (FM27b) is **not**
  PSD on the whole lower branch.  At `(m,N,t) = (200,1842,4)`,
  `Q_2(0,1) < 0`, while the actual target is positive.  Exact checks of 33
  lower-branch and 110 upper-branch points are all positive.
- **FM-T1-2L6 (saturn).**  The correction becomes a pair of binomial tails:
  `M_2 = -sum_(l,w) eta_l det(w) [ sum_(u>=1) b_a(alpha+u) C_(a,t)(beta-u) + sum_(u>=1) b_a(beta+t+u) C_(a,t)(alpha-t-u) ]`.
  Definitions (added after FM-CHK9):
  - `b_a(x) = binom(a,x)` for integers `0 <= x <= a`, and 0 otherwise;
  - `C_(a,t)(x) = sum_(h=1)^t b_a(x+h)`;
  - `l = 0,1,2` with `y_0 = (2,1)`, `y_1 = (3,2)`, `y_2 = (4,1)` and
    `eta = (5,-3,1)`;
  - `w` runs over `W(C_2)` (signed permutations), with `(A,C) = w y_l`;
  - `alpha = (a+s+t+3-A-C)/2` and `beta = (a-s-t-1+A-C)/2`;
  - `s >= t`, `s+t+a` even.

  Checks: `M_2(1,0,a) = 0`, `M_2(1,1,0) = -3`, `M_2(2,1,1) = -2`.
  - The full Vandermonde convolutions cancel over `W(C_2)`.
  - `t = 6` values: `L(7,6,1) = 2`, `L(8,6,2) = 2`, `L(9,6,3) = 3`.
  - No bound yet.
- **FM-L123f (neptune): proved** the `m = r = 2` inequality
  `3m20' <= m22' + 2m11' + 6m00'` for `G = M_(2k) M_(2l)` (two even minus
  generators, all `k`, `l`).  The margins are:

  | case | margin |
  |---|---|
  | `a = b = 0` | 6 |
  | `a = b >= 1` | 4 |
  | `|a-b| = 1` | 2 |
  | `|a-b| = 2` | 1 |
  | otherwise | 0 |

  Here `a = k-1` and `b = l-1`, via SO(5) harmonics `V'_(b,b) = [b,0]`.
  - Exact screens are clean: 792 multisets of up to five even factors, and
    a further 2,782 mixed multisets.
  - Raw Brauer diagrams are linearly dependent in dimension 5 (the 6x6
    determinant relation), so a naive diagram injection cannot work.
- **Turán generating function (main agent).**  With
  `c_n = [z^n](1+z)^a(1-z)^t`:

  ```text
  D_k = sum_(n even) Cat(n/2) c_n binom(N-n, k-n/2),
  G(x) = sum_k D_k x^k = (1+x)^N sum_j Cat(j) c_(2j) (x/(1+x)^2)^j,
  G(x) = < F(sqrt(x) lambda) F(sqrt(x)/lambda), 1 >_(SU(2)).
  ```

  - The identity holds exactly in 6,000 cases.
  - The right-half monotonicity `D_k >= D_(k+1)` is exactly unimodality of
    the symmetric coefficient sequence of `G`.
  - It would follow from `c_(2j) >= 0` (`gamma`-positivity), which fails in
    general, or from real-rootedness of `G`, which is being tested.

### FM-CHK-B1 (luna_max_venus, second independent checker)

- **FM29: ACCEPTED.**  The checker gave a symbolic re-derivation of the
  companion recurrence through a one-particle moment recursion, gave
  explicit `P_3`, `Q_2`, `Q_3` with all coefficients positive (minima 27,
  45, 135), and verified the Cauchy--Schwarz argument.
- **FM32: ACCEPTED.**  An exact general formula for `E[f]` through the
  sequence `p_k` gives certificates for `E[S^4 G]` and `E[S^4 G H_3]`
  (smallest coefficients 10368 and 103680).  `E[Q_5]` and `E[S Q_5]` are
  proved by grouping; both vanish exactly at `(e,r) = (1,0), (2,0)`.
- **FM34: ACCEPTED.**  `P_3 H_3 = V'_22 + V'_11 + V'_00`, with dimensions
  `14 + 5 + 1 = 20`.  The FM5 absorption needs the case split `q = p` (no
  residual) versus `q > p`.
- **FM35: ACCEPTED.**  The checker gave an independent symbolic proof
  through the closed-form walk count
  `P_a(d_1,d_2) = a!^2/((M+d_1)!(M+d_2)!(M+d_1+d_2)! M!)`, including all
  boundary cases.
- **FM30 extension: sign argument ACCEPTED.**  The checker flagged the
  identity for `M` as needing a symbolic proof.  That proof is already
  recorded in the FM-CHK7 repairs (main agent, sympy, identically in
  `(N,d,e)`), so FM30 and its extension are complete.

### FM-T2-PQ7 (luna_max_jupiter): two plus labels, nearly complete; one certificate defect

**Verified by the main agent.**

- *Exact screen.*  All `1,045,396` negative-window triples with
  `4 <= q < 50`, `p < 3q`, `h < 2A` were re-run with exact rational
  arithmetic.  The maximum of `R B + U B^2` is
  `1916517521177/4962905706564 ~ 0.386`, at `(4,5,5)`.  With FM-T2-PQ6
  (`q = 2, 3`), the family is proved for all `q <= 49`.
- *Lemmas for `q >= 50`.*
  - `Delta(2q,q) = -2q^4 + 24q^3 + 54q^2 + 36q + 9 < 0`, and the
    `d Delta/dp` bound (grid check), so `p < 2q`.
  - `|C/q^4 - (x^2-1)^2| <= 18.5 z <= 19 z`.
  - The bound `-Pi/q^4 <= 8t(1-t) + 19z` (completing the square).
  - The small-suffix case `h <= q`: `60 q^3 e^(-3q/8) < 1`.
- *Long-suffix bounds (2), (3):* no violation on 373,354 random admissible
  samples with `q >= 50`.

**Defect.**  The rational interval certificate for the long suffix
(`h >= q`) bounds `(x+1+2z)^4` by `(51/25)^4 = 2.04^4`.  But `x = p/q` ranges
up to 2, so the correct constant is `(76/25)^4 = 3.04^4`.  With the corrected
constant the top bin (`1/10 <= t <= 1.05`) gives `1.52 > 24/25`, and the
certificate fails as written.  The claim itself looks true: the actual
maximum of `R B + U B^2` over 300,000 random long-suffix samples with
`q >= 50` is `0.339`.

**Status.**  The two-plus-label family is proved except for `q >= 50` with
`h >= q`, which needs a repaired certificate (FM-T2-PQ8).

**Kill K21 (Turán generating function not real-rooted; main agent, Sturm counts).**
`G(x) = sum_k D_k x^k` is not real-rooted in 192 of 222 cases with
`a < 25`, `t <= 8`.  The smallest case is `a = t = 1`: `c = (1,0,-1)` gives
`G = 1 + x + x^2`.  Unimodality of `(D_k)` therefore cannot come from
real-rootedness plus Newton.  It must come from the `gamma`-type /
Catalan structure, from an injection, or from the explicit residual
analysis.

### One-label residual: status after FM-T1-AB10 (mars) and FM-T1-UNI (venus)

- **mars:** gave an analytic proof of the Catalan identity, and the exact
  bracket form
  `D_k - D_(k+1) = sum_j Cat(j) c_(2j) b(N-2j, m)`, where
  `b(R,m) = 2(m+1)/(R+m+2) binom(R,(R+m)/2)` is the ballot number
  `int U_m(x) x^R dsigma`.  The bracket sequence is log-convex at some
  indices, so a log-concavity pairing fails.
- **venus:** real-rootedness fails at `(a,t) = (1,1)` (kill K21, found
  independently).  Coefficientwise `gamma`-positivity fails there too
  (`c_2 = -1`).  `D` is not log-concave, e.g. `(a,t) = (3,11)`.  Exact
  unimodality holds on all 496 pairs `0 <= a <= t <= 30`.
- **main agent:** summing the ballot form back gives
  `D_k - D_(k+1) = int int U_m(x) (x+y)^a (x-y)^t dsigma dsigma`, the
  original single-label integral.  So the Catalan identity is a
  reformulation and adds no new positivity.

**Residual.**  The branches `R_-` and `R_+` (labels `>= 6`, `t >= 4`)
remain open.  They are recorded here with all tools, and the line is
paused.  [Closed on 2026-09-30 by Theorem OL, item (39).]

### FM-CHK9 (luna_max_mercury) and FM-L123g (luna_max_neptune)

**FM-CHK9 (mercury).**

- *ACCEPTED:*
  - the FM40 branching identity;
  - neptune's two-even-factor theorem, with a direct `B_2` derivation of
    `[a,0] (x) [b,0] = sum_(j,k) [a+b-2j-k, k]`;
  - jupiter's `q = 2, 3` windows (37 cases), with 10 values independently
    recomputed;
  - the Catalan identity, with a proof.
- *Repairs, now applied:* the FM40 parity clause, the definitions for the
  two-tail formula, and the ledger scope.

**FM-L123g (neptune).**  The half-space `Delta(W) >= 0` on arbitrary
multiplicity vectors is not preserved by `[1,0]`.  Example:
`Delta([2,1]) = 0` but `Delta([2,1] (x) [1,0]) = -2`.  The word-level
inequality is not contradicted, since `[2,1]` is not a word module.  The
first open word case is three even factors `[b_1,0] (x) [b_2,0] (x) [b_3,0]`.

### Corollary FM41 (two symmetric powers, smaller one at most 7; FM-T1-2L7, luna_max_saturn; accepted)

For `t = 6, 7`, with `n = s-t`, `a = n+2v`:

```text
L(s,t,a) = binom(n+2v,n+v)^2 (s+1)(s+2)(s+3) F_t(s,v) Pi_t(n,v) / D_t(s,v),
F_6 = (s+2v-5)(s+2v-4),   F_7 = (s+2v-6)^2 (s+2v-5)(s+2v-4),
Pi_t = sum_(j=0)^14 A_(t,j)(n) v^j.
```

Coefficient tables (inserted after FM-CHK10).  Here
`P(c_0,...,c_d) = sum c_i n^i`, with coefficients in ascending powers.

```text
t=6 (n = s-6):
A_6,0  = (n+1)(n+2)(n+3)^2(n+4)^2(n+5)^2(n+6)^2(n+7)(n+8)(n+9)(n+10)(n^2-3n+6)
A_6,1  = 2(n+3)(n+4)(n+5)(n+6) P(42575040,96146304,80907816,34490004,10637178,3905387,1478946,388869,63102,6069,318,7)
A_6,2  = 4(n+3) P(5203091520,12942002352,13520571544,7999290884,3158730314,1004135756,308577253,89770777,20996567,3529447,402361,29401,1241,23)
A_6,3  = 8 P(9643348080,23660009424,25311804720,15739314810,6501823266,2022370649,552406008,143648452,32575362,5631108,677484,52726,2376,47)
A_6,4  = 2 P(29681520480,64864655736,62356073124,34763522606,12719704579,3448679916,815618904,183649854,35359305,4963132,452162,23604,534)
A_6,5  = 12 P(2598530760,5208553056,4537687044,2240890949,705851652,161919914,32984680,6507647,1052344,114096,7048,186)
A_6,6  = 4 P(3208721040,5851087998,4481667597,1874058292,485258120,92092608,16525818,2887170,374367,28040,886)
A_6,7  = 24 P(187071030,290218584,184516872,62489398,13026708,2086460,338376,49847,4614,181)
A_6,8  = 6 P(198903480,256749276,131476058,35988094,6160155,837544,116598,13222,685)
A_6,9  = 12 P(18779580,19553628,8344352,1906673,253964,24890,3080,245)
A_6,10 = 12 P(2764380,2265542,837193,143518,8765,504,126)
A_6,11 = 504 P(4410,5548,1568,109,-10,1)
A_6,12 = 84 P(7740,3058,491,-46,1)
A_6,13 = -1008(n^2-19n-15),   A_6,14 = 5040
t=7 (n = s-7):
A_7,0  = (n+2)(n+3)(n+4)^2(n+5)^2(n+6)^2(n+7)^2(n+8)(n+9)(n+10)(n+11)(n^2-3n+6)
A_7,1  = (n+4)(n+5)(n+6)(n+7) P(500083200,816663744,519495864,170869788,46683682,17484515,5972107,1320170,178588,14411,639,12)
A_7,2  = P(444072153600,1031536195776,1045050403728,612945866736,241097310836,73973451540,20793147450,5566614039,1272340089,222441051,28001461,2434869,138447,4629,69)
A_7,3  = 2 P(187139320704,384520036176,342701018904,177249811476,61470187566,16573326788,4098471489,967609969,192742713,28563339,2920407,192991,7401,125)
A_7,4  = P(213542369856,392103163584,314172267348,145894089948,44989136057,10611471696,2271103518,462565704,78111321,9427452,734272,32784,636)
A_7,5  = 6 P(14109635952,22966671132,16570921540,6892460691,1864752381,379186077,70590181,12561290,1775700,165890,8806,200)
A_7,6  = 2 P(12454768728,19133033922,12422834487,4452982881,1000180111,168546615,27544815,4318188,489693,31734,866)
A_7,7  = 12 P(502172766,680756178,381479546,112895968,20442248,2916856,439385,58738,4767,162)
A_7,8  = 3 P(485206104,522123156,226819174,52386516,7609405,941904,125038,12672,567)
A_7,9  = 6 P(27117468,28719168,10843776,2122999,253149,24943,2847,190)
A_7,10 = 6 P(5699964,3131370,902963,132843,8500,567,93)
A_7,11 = 36 P(-4530,30702,8378,611,-40,5)
A_7,12 = 6 P(67332,16758,2551,-210,5)
A_7,13 = -360(n^2-17n+9),   A_7,14 = 1800
```

The `n = 0` split for `t = 7`: the last five terms are `v^10` times
`34199784 - 163080v + 403992v^2 - 3240v^3 + 1800v^4`, which equals
`(34199784 - 163080v + 201996v^2) + v^2 (201996 - 3240v + 1800v^2)`.  The
discriminants are `-27606283189056` and `-1443873600`.

**Checks (main agent).**

- *Identities:* both closed forms are exact symbolic identities in `(s,v)`,
  obtained through the constituent sum of the 24-term formula with
  binomial shift ratios (93 s and 150 s in sympy).
- *`t = 6` certificate:* Sturm counts show `A_(6,j) > 0` on `[0, inf)` for
  `j <= 12`.  The final quadratic `A_(6,12) + A_(6,13) v + A_(6,14) v^2`
  has discriminant with no root in `[0, inf)` (value `-12878611200` at 0).
- *`t = 7` certificate:*
  - `A_(7,j) > 0` for `j <= 10` and for `j = 12`.
  - `A_(7,11)` has its only real root at `n ~ 0.142`, so it is positive for
    integers `n >= 1`.
  - At `n = 0`, the packet's quartic split works: both quadratics have
    negative discriminants.
  - The final-quadratic discriminant has no root in `[0, inf)`.

**Corollary FM41.**  Q3 holds for minus labels `{s+1, t+1, 1, 1}` with any
number of plus 1s whenever `min(s,t) <= 7`.  Open: `s > t >= 8` (interior).

### Theorem FM42 (two arbitrary plus labels with a fundamental suffix; FM-T2-PQ2..PQ8, luna_max_jupiter; accepted)

For all `p, q >= 1` and `a >= 0`:

```text
< Shat_p Shat_q (C^4)^(x)a, 1 >_(Sp(4)) >= 0.
```

Equivalently, Q3 holds for every word `[1^-, 1^-, p^+, q^+, 1^+ x a]`.

**Notation for the certificate** (added after FM-CHK10).  Put
`h = (a+2-p-q)/2` and

```text
N = p+q+2h,   S = p+q+h+2,
D = 3N^2 + (p-q)^4 - 4(p-q)^2,   E = 3N^2 + (p+q+2)^4 - 4(p+q+2)^2,
A = p^2+q^2-p-q-3,   C = (p+q+1)(p+q+2)((p-q)^2-p-q),   Pi = 12h^2 - 4Ah + C,
R = N(p+1)(q+1)(-Pi)/(h S D),   U = E(h+p+1)(h+q+1)/(h S D),
rho = prod_(i=0)^q (h+i)/(p+q+h+1-i) <= B = ((2h+q)/(2h+2p+q+2))^(q+1).
```

With these, `T = alpha (1 - R rho - U rho^2)` and `alpha > 0` when
`Pi < 0`.

*Small-suffix bounds (`q >= 50`, `h <= q`).*

- `A < 5q^2`, and when `C < 0`, `|C| < 30q^3`; hence
  `-Pi <= 20q^2 h + 30q^3`.
- `E <= 87q^4`, `h+p+1 <= 3q`, `h+q+1 <= 2.02q`, `S >= 2q`, `D >= 2N^2`,
  `N >= 2q`.  Substitution gives `R <= 26q^3` and `U <= 33q^3`.

**Proof assembly.**

1. *Parity.*  `T = 0` when `a+p+q` is odd.
2. *Structural cases.*  The telescoping
   `T = Delta_d - Delta_(e+1) + M` holds, with `M` in closed form (FM30,
   proved as a symbolic identity).  This gives:
   - `p+q >= a+2` (boundary and outer support);
   - `p = q`;
   - `p` or `q` equal to 1 (FM24);
   - `a = 0`;
   - `Pi >= 0` (`M >= 0` there, and FM24 positivity).
3. *`q = 2, 3` with `Pi < 0`.*  The 37 window cases were computed exactly
   (FM-T2-PQ6, checked by FM-CHK9).
4. *`4 <= q <= 49` with `Pi < 0`.*  The Jensen certificate
   `R B + U B^2 < 1` holds on all 1,045,396 triples, maximum 0.386 (main
   agent re-ran the exact screen).
5. *`q >= 50`, `h <= q`.*  `R B + U B^2 <= 60 q^3 e^(-3q/8) < 1`.
6. *`q >= 50`, `h >= q`.*  A rational interval certificate
   (FM-T2-PQ8):
   - `t <= 1/50`: monotonicity in `t`, value `0.00796`;
   - `1/50 <= t <= 1/10`: 80 bins, `x = 1` worst case since the
     `x`-factors decrease in `x` for `t + eps < 1`, maximum `0.268`;
   - `1/10 <= t <= 21/20`: 38,000 two-dimensional `(t,x)` rectangles, upper
     `x`-endpoint in polynomial factors and lower endpoint in exponentials,
     maximum `0.8844 < 9/10`;
   - `e^(-y) <= (1+y/128)^(-128)` keeps everything rational.

**Checks.**

- The main agent re-ran the verifier exactly (all three assertions pass) and
  checked the monotonicity-in-`x` and bin-endpoint arguments.
- The earlier constant defect (2.04 versus 3.04) no longer matters.  The
  2.04 bound is now used only where `x = 1` is provably worst.
- The independent check of the `q >= 50` lemmas is pending (FM-CHK10).

This is a complete infinite family with **two arbitrary labels**, and the
second complete family of arbitrary labels after FM24.

### Corollary FM44 (three even minus labels at `m = r = 2`; FM-L123i, luna_max_neptune; accepted pending check)

Put `Q = [2,0] - 3[1,1] + 2[1,0] + 6[0,0]` (SO(5)), and write `R(a,b)` for
the multiplicity-free support of `[a,0] (x) [b,0]`.

**Claim.**  `Delta(H_a H_b H_c) = <H_a H_b H_c, Q> >= 0` for all `a, b, c`.

**Proof** (packet).  By Frobenius, `Delta` becomes an indicator sum over
`R(a,b)`, with exact case tables for `c = 0, 1`.  For `c >= 2` it is formula
(2) of the packet.  Each negative level-1 indicator forces enough positive
level-0 indicators:

- `eta_-` forces `4`;
- `eta_0` forces `2`;
- `eta_+` forces `5`.

Parity separates `eta_0` from `eta_+-`.

**Words covered.**  Minus labels `{2a+2, 2b+2, 2c+2, ...}` at `m = r = 2`.

**Open.**  Four or more factors.  For a prefix `P`, the exact local-band
inequality (3) is needed.  Constituentwise induction fails (the `[2,1]`
example).

### FM-CHK10 (mercury)

- The FM40 two-tail formula is **ACCEPTED**.  Its proof is the Vandermonde
  cancellation over `W(C_2)`, compared in 405 cases.
- The `q >= 50` lemmas of FM-T2-PQ7 are **ACCEPTED**: `p < 2q` by a
  derivative proof, the `C`-bound `18.5z`, the `-Pi` bound, and the
  small-suffix `B` bound.  The `R`, `U` definitions it asked for are now
  inserted above.
- The `t = 6, 7` tables, which were missing when the check ran, are now
  inserted.
- The long-suffix defect was confirmed independently.  It is repaired in
  FM-T2-PQ8 (FM42).

### Corollary FM45 (two odd minus labels at `m = r = 2`; FM-L123h, luna_max_venus; accepted)

Let `I_kl = int H_(2k+1) H_(2l+1) u^4 v^4 dsigma dsigma`.  Put
`F_k = (x^2-y^2)^2 H_(2k+1)`.  Then `I_kl = <F_k, F_l>` in the orthonormal
basis `U_i(x) U_j(y)`, and both coefficient vectors are explicit.

The values:

| case | `I_kl` |
|---|---|
| `k, l >= 1` | 30, 18, 10, 2, 0 for `|k-l| = 0, 1, 2, 3, >= 4` |
| `k = l = 1` | 38 (exception) |
| `k = 0` | `I_00 = 12`, `I_01 = 16`, `I_02 = 6`, `I_03 = 2`, else 0 |

All are `>= 0`.  The main agent recomputed the full table exactly.

**Words covered.**  `[2k+1^-, 2l+1^-, 1^-, 1^-, 1^+ x 4]`, i.e. the
`a = 4`, even-`(s,t)` slice of the two-large-parts family, for all labels.

### FM-T2-PQR (luna_max_mars): three plus labels, partial

**Exact decomposition.**

```text
T(p,q,s;a) = sum_l A_l tau_l(a) + sum_(j in C(p,q)) M(j,s;a)
           + sum_(j in C(p,s)) M(j,q;a) + sum_(j in C(q,s)) M(j,p;a),
```

where `tau_l` is the FM24 one-label value and `M` is the FM30 cross term
(with a moment formula).  Grouped form:

```text
T = sum_(j in C(p,q)) T(j,s;a) + (the other two correction sums).
```

**Proved.**

- Whenever the corrections vanish by the support rule.
- The grouped condition `(6)` via FM30 and FM-T2-PQ6.
- One label equal to 1, which reduces to FM42.
- `a <= 2`, by the seven-factor theorem.
- Labels in `{1,2}`.

**Checks.**  `T(1,1,1;0) = 0`, `T(1,1,2;0) = 1`, `T(3,3,2;0) = 6`,
`T(1,1,1;1) = 3`, all matching the main agent.

With FM42 now complete, every grouped term `T(j,s;a) >= 0`; the correction
sums are what remain.

### Corollary FM43 (two symmetric powers, smaller one at most 9; FM-T1-2L8, luna_max_saturn; accepted)

For `t = 8, 9`, with `n = s-t` and `a = n+2v`:

```text
L(s,t,a) = binom(n+2v,n+v)^2 (s+1)(s+2)(s+3) F_t Pi_t / D_t,
F_8 = (a+1)(a+2),   F_9 = (a+1)^2 (a+2)(a+3),
Pi_t = sum_(j=0)^18 A_(t,j)(n) v^j,
```

with `D_t` as in FM41.  Tables, where `P(c_0,...)` is ascending in `n`:

```text
t = 8

A_8,0 = (n+1)(n+2)(n+3)^2(n+4)^2(n+5)^2(n+6)^2(n+7)^2(n+8)^2
        (n+9)(n+10)(n+11)(n+12)(n^2-3n+6)
A_8,1 = 2(n+3)(n+4)(n+5)(n+6)(n+7)(n+8)
        P(6250279680,15146260704,14236651800,7131840276,2510259942,
          896031797,342708629,105509586,22504611,3213030,301611,17878,607,9)
A_8,2 = (n+3)P(186140420720640,548455281596928,710149171034880,
          544244943756672,284150047037904,113370979615968,38918368007824,
          12411105478644,3596573229371,878461217611,170492445705,
          25462296317,2867844971,238483119,14177491,569627,13854,154)
A_8,3 = 2P(375373017704448,1108296861716736,1476626629606752,
          1184271710014048,648139296826944,264621929399660,88173776645646,
          26074204221205,7055098697746,1678224270170,329533656580,
          50925872012,6001199112,524634410,32842934,1389807,35582,416)
A_8,4 = P(655944971584896,1781490276819360,2192433630663888,
          1621089546453280,812769319676792,301187328559160,90185748508568,
          23799814798062,5726942501707,1204717706604,206666153748,
          27364527382,2686627900,187838400,8815832,248600,3181)
A_8,5 = 2P(206642024831520,523142182925760,595001438706216,
          402391550266464,182587584271378,60679625592648,16210562715848,
          3816782354426,818520196470,151802012344,22477697246,
          2492621253,196438240,10325062,323658,4571)
A_8,6 = 2P(103192513698240,240023090871096,247853825300136,
          150662048954088,60986759516543,17992708012244,4266396881601,
          895240233472,170692516675,27634243662,3458949907,309933790,
          18449124,650048,10238)
A_8,7 = 4P(20902989591288,43507088259972,40046854930104,
          21639634647273,7758103979396,2017818790849,421409997514,
          78292183781,13170679662,1830895398,187515556,12799937,513848,9142)
A_8,8 = P(26383558068792,49105241363730,40272033105423,
          19336462949538,6100965968947,1378374297036,249469491176,
          40777732386,6020784744,701560356,55618688,2586198,52834)
A_8,9 = 2P(3282608199330,5449905930990,4011731361282,1702089398210,
          462744465514,88234109810,13650753462,1982725461,255140430,
          23676562,1299650,31141)
A_8,10 = 2P(664021608390,999889896732,651000162354,236643112204,
           53285372161,8347751040,1119555171,147047886,15745065,1038818,30041)
A_8,11 = 132P(1662806844,2305866284,1284774668,383416645,69380826,
           8993733,1078570,123320,9828,358)
A_8,12 = 66P(567620208,575786390,258776965,61513002,9011652,
           973634,99278,9122,453)
A_8,13 = 132P(21580230,26383850,9503314,1882927,215948,16289,1344,112)
A_8,14 = 132P(3858810,1765810,579141,86270,5422,60,41)
A_8,15 = 264P(-13854,59068,13226,1233,-80,5)
A_8,16 = 33P(111048,21490,4783,-310,5)
A_8,17 = -1980(n^2-25n+17)
A_8,18 = 9900

t = 9

A_9,0 = (n+2)(n+3)(n+4)^2(n+5)^2(n+6)^2(n+7)^2(n+8)^2(n+9)^2
        (n+10)(n+11)(n+12)(n+13)(n^2-3n+6)
A_9,1 = (n+4)(n+5)(n+6)(n+7)(n+8)(n+9)
        P(88333632000,158822354304,115262715168,45605264472,
          14048920956,5013505842,1794202011,475981340,85684605,
          10348818,828021,42168,1239,16)
A_9,2 = P(6670275315148800,18106834387548672,22075286669359488,
          16136618995623456,8093681318686704,3098292449366736,
          1007643964900192,300930866289174,82327045581507,
          19433890472889,3742551364409,569051379210,67054479128,
          6027863712,404896788,19660560,651661,13191,123)
A_9,3 = P(6244758017519616,15486959430660096,17299375657004160,
          11645222362677376,5385584726557552,1889558130947124,
          556206055876656,148734944803055,36355017533715,7669562288061,
          1313521830335,175511505409,17831482269,1344255553,
          72687693,2662916,59172,602)
A_9,4 = P(4025380715796096,9205929207666336,9542004587136704,
          5960042110910736,2543474889068700,814523002475840,
          215833434592192,51460243083869,11188561210318,2093468440661,
          314815993151,36264392595,3089811567,187514439,7650229,187908,2099)
A_9,5 = P(1897708847165376,4012149575153312,3856287017119072,
          2224097636672044,868414103987208,251674012064905,
          59843275922081,12770682852521,2480263962465,410321508828,
          53471918742,5183199649,356298861,16339817,447507,5532)
A_9,6 = P(696105392743680,1387916738373632,1231385266503980,
          645411802548068,226339536031441,58449786188983,
          12369553659679,2359998593647,408674695436,59161486056,
          6535524732,513642027,26772451,826371,11417)
A_9,7 = 2P(106725147085872,192134650968196,153941869011750,
          72370374912952,22603535186338,5172246385714,971134287548,
          165324522375,25350408783,3153131367,285268830,17115251,603351,9433)
A_9,8 = P(55174467453480,88498768689422,62573742396361,
          25969295103238,7128225647984,1422922753002,233505345426,
          35177934432,4726999839,491379606,34399332,1408044,25314)
A_9,9 = P(10103910079788,15104855388308,9784720935800,3635483042492,
          870182867364,148378743165,21031281015,2821258158,328877898,
          27163393,1317151,27808)
A_9,10 = P(1972144573380,2364237988112,1307239302548,413846136616,
           82272716243,11560449703,1425615425,173084230,16769701,983851,25091)
A_9,11 = 22P(7826024472,11872950220,6105489206,1634296520,262786446,
           30597073,3425152,365840,26268,843)
A_9,12 = 11P(3830536992,2812956500,1068494690,220429016,28166875,
           2772892,274274,23272,1009)
A_9,13 = 44P(14594454,45780080,15601185,2737340,280541,21063,1708,119)
A_9,14 = 44P(10296810,3200696,827170,105073,6512,147,42)
A_9,15 = 88P(-140238,66274,13948,1311,-70,5)
A_9,16 = 11P(184392,23762,4891,-290,5)
A_9,17 = -660(n^2-23n+49)
A_9,18 = 3300
```

**Checks (main agent).**

- *Identities:* both are exact symbolic identities in `(s,v)` (171 s and
  250 s).
- *Coefficient signs (Sturm counts):* `A_(t,j) > 0` on `[0, inf)` for
  `j <= 14` and `j = 16`.
- *`A_(8,15)`:* its only real root is at `n ~ 0.223`, so it is positive for
  `n >= 1`.  At `n = 0` the quartic split has discriminants
  `-3719829429258624` and `-71425767600`.
- *`A_(9,15)`:* its only real root is at `n ~ 1.546`, so it is positive for
  `n >= 2`.  At `n = 0` and `n = 1` the splits have discriminants
  `(-1685593710244224, -12340983600)` and
  `(-2946464029238400, -15128823600)`.
- *Final quadratic:* `A_(t,16) + A_(t,17) v + A_(t,18) v^2` has discriminant
  with no root in `[0, inf)`; its value at 0 is `-143984530800` (`t = 8`)
  and `-25727842800` (`t = 9`).

**Corollary FM43.**  Q3 holds for minus labels `{s+1, t+1, 1, 1}` with any
number of plus 1s whenever `min(s,t) <= 9`.

**Patterns for `t = 1..9`** (saturn):

- `deg_v Pi_t` is `2t+2` for even `t` and `2t` for odd `t`;
- the leading coefficient is `20(k+1)(2k+1)(2k+3)` for `t = 2k` and
  `10(k+1)(k+2)(2k+3)` for `t = 2k+1`;
- `F_t` is `(a+1)(a+2)` for even `t` and `(a+1)^2(a+2)(a+3)` for odd `t`
  (`t >= 4`).

The case of all `t` is FM-T1-2L9.

**One-label residual, final status for this session.**  Three independent
attempts (mars FM-T1-AB8..AB10, venus FM-T1-UNI, jupiter FM-T1-RES) did
not close the residual branches `R_-` and `R_+`.  Labels `>= 6`, both
`a, t >= 4`, and `F_m > 0` remain open.  The exact screens are all positive
(96,929 cases plus the branch points).  The line is paused.

### FM-T1-2L9 (luna_max_saturn): exact all-`t` formula for two symmetric powers

Put `C_* = binom(n+2v, n+v)` and `R_delta = binom(a, n+v+delta)/C_*` (exact
shift ratios).  With `g_q = sum_r (-1)^r binom(3,r) R_(q-r)` and
`d_q = g_q^2 - g_(q-1) g_(q+1)`:

```text
L(s,t,a)/C_*^2 = d_2 - d_(t+3) - sum_(l,w) eta_l det(w) (T1_(l,w) + T2_(l,w)),
```

where `T1` and `T2` are explicit finite double sums of products of `R`s (in
the packet).  This holds for every `t`, and the one-label axis part
telescopes exactly for all `t`.  Positivity for all `t` is the remaining
comparison (3) (tails versus `d_2 - d_(t+3)`).  Open for `t >= 10`.

### Sign map of `Lambda(p,q;a) = W_2(V_(p,q) (x) (C^4)^(x)a)` (main agent, exact)

Negative values occur only for constituents with `q >= 1`, and only in a
band.  The negative entries by `a`:

| `a` | negative `Lambda` |
|---|---|
| 0 | `(1,1): -3` |
| 1 | `(2,1): -2` |
| 2 | `(2,2): -2`, `(3,1): -1` |
| 3 | `(3,2): -3` |
| 4 | `(2,2): -3`, `(3,1): -1`, `(3,3): -3`, `(4,2): -3` |
| 6 | `(3,3): -16`, `(4,2): -18`, `(4,4): -6`, `(5,3): -8` |
| 8 | `(3,3): -72`, `(4,2): -81`, `(4,4): -62`, `(5,3): -91`, `(5,5): -14`, `(6,2): -10`, `(6,4): -22`, `(7,3): -5` |

The pure symmetric constituents `(p,0)` are never negative (FM19).  Every
H-only word at `r = 2` is a nonnegative combination
`sum_lambda m_M(lambda) Lambda(lambda; a)` with `M` a product of symmetric
powers.  So the core is a multiplicity-domination statement for tensor
products of symmetric powers.  It generalises T1's `m_11 <= m_20`, which
the repository proved through the 3-nesting-free Pfaffian model.  Assigned
as FM-T1-DOM.


### Appendix to FM42 (long-suffix certificate, verbatim from FM-T2-PQ8; inserted after FM-CHK11)

Notation: `x = p/q`, `t = h/q^2`, `z = 1/q`, `P = -Pi/q^4`,
`eps = (5/2) z + z^2 <= 63/1250`, and `k_0 = 2601/14995`.

#### A.2 Bounds in the negative-`Pi` region, `q >= 50`
Bounds in the negative-\(\Pi\) region for \(q\ge50\)

Write
\[
A=p^2+q^2-p-q-3,\quad
C=(p+q+1)(p+q+2)((p-q)^2-p-q),\quad
\Pi=12h^2-4Ah+C.
\]
The discriminant in \(h\) is \(16\Delta\), where \(\Delta=A^2-3C\). Expanding,
\[
\begin{aligned}
\Delta={}&-2p^4-8p^3+8p^2q^2+16p^2q-2p^2+16pq^2+32pq+12p\\
&-2q^4-8q^3-2q^2+12q+9.
\end{aligned}
\]
For \(p\ge2q\), \(q\le p/2\), so
\[
\partial_p\Delta\le-4p^3-4p^2+12p+12<0
\]
because \(p\ge100\). Also
\[
\Delta(2q)=-2q^4+24q^3+54q^2+36q+9<0\quad(q\ge50).
\]
Since \(\Pi<0\) requires \(\Delta>0\), this proves \(p<2q\).

Put \(x=p/q\), \(t=h/q^2\), \(z=1/q\), and \(P=-\Pi/q^4>0\). Then \(1<x<2\), \(z\le1/50\), and
\[
\frac A{q^2}=x^2+1-(x+1)z-3z^2,
\]
\[
\frac C{q^4}=(x+1+z)(x+1+2z)((x-1)^2-(x+1)z).
\]
The coefficients of \(z,z^2,z^3\) in the latter expression minus \((x^2-1)^2\) have absolute values at most \(18,25,6\), respectively. Hence
\[
\left|\frac C{q^4}-(x^2-1)^2\right|
 \le18z+25z^2+6z^3\le19z.
\]
Completing a square now gives
\[
0<P\le -12t^2+4(x^2+1)t-(x^2-1)^2+19z
 \le8t(1-t)+19z. \tag{1}
\]
If \(t\ge1\), (1) implies \(t-1<19/(8q)\le19/400\). Thus \(t<21/20\). Also \(h\ge q\) implies \(t\ge z\).

For the Jensen bound, the function
\[
v\longmapsto\log\frac{h+v}{h+p+1+v}
\]
is concave. Applying Jensen to \(v=0,\ldots,q\) gives \(\rho\le B\). Further,
\[
B\le \exp\!\left(-\frac{x}{t+\varepsilon}\right),
\qquad \varepsilon=\frac52z+z^2\le\frac{63}{1250}. \tag{2}
\]
For completeness, the exponent in \(B\le e^{-\lambda}\) is
\[
\lambda=\frac{2(1+z)(x+z)}{2t+2xz+z+2z^2}.
\]
After clearing positive denominators in \(\lambda\ge x/(t+5z/2+z^2)\), the difference has the sign of
\[
2tz(x+1+z)+z\bigl(2x(2-x)+(5x+5)z+(2x+7)z^2+2z^3\bigr),
\]
which is nonnegative for \(1\le x\le2\).

Using \(N\ge2q^2(t+z)\), \(D\ge3N^2-4\ge2.999N^2\), and \(N/S\le2\), direct substitution gives
\[
R\le k_0\,\frac{xP}{t(t+z)^2},\qquad
U\le\frac{1+2z/t}{2.999}
 \left(3+\frac{(x+1+2z)^4}{4(t+z)^2}\right),
\quad k_0=\frac{2601}{14995}. \tag{3}
\]
Here \(D\ge2.999N^2\) follows from \(N\ge2h\ge100\). The factor \((x+1+2z)^4\) is retained in the two-dimensional certificate below.

#### A.4 Exact Fraction certificate for `h >= q`, `q >= 50`
Exact Fraction certificate for \(h\ge q\), \(q\ge50\)

When \(z\le t\le1/50\), (1) implies \(P/[t(t+z)^2]\le13/t^2\). The \(x\)-dependent factors in (3) times their exponentials are decreasing for \(x\ge1\) on this interval. Since \(t+\varepsilon\le3.52t\), the bound is at most
\[
k_0\frac{13}{t^2}e^{-1/(3.52t)}
+\frac3{2.999}\left(3+\frac{2.04^4}{4t^2}\right)e^{-2/(3.52t)}.
\]
Each term increases for \(t\le1/50\), and its endpoint value is \(<1/100\).

For \(1/50\le t\le1/10\), the \(x\)-dependent factors are again decreasing for \(x\ge1\). I used 80 rational bins of width \(1/1000\). For \(a\le t\le b\), the bin bound is
\[
\begin{aligned}
&k_0\left(\frac{8(1-a)}{a^2}
 +\frac{19/50}{a(a+1/50)^2}\right)\mathcal E\!\left(\frac1{b+\varepsilon_0}\right)\\
&\quad+\frac{1+1/(25a)}{2.999}
 \left(3+\frac{(51/25)^4}{4a^2}\right)
 \mathcal E\!\left(\frac2{b+\varepsilon_0}\right),
\end{aligned}
\]
where \(\varepsilon_0=63/1250\) and \(\mathcal E(y)=(1+y/128)^{-128}\).

For \(1/10\le t\le21/20\), I used 190 bins in \(t\) and 200 bins in \(x\), both of width \(1/200\). On a rectangle \(a\le t\le b,\ c\le x\le d\), put
\[
G(a)=\begin{cases}8(1-a)/a^2,&a<1,\\0,&a\ge1,\end{cases}
\qquad C(a)=\frac{19/50}{a(a+1/50)^2}.
\]
The exact rational upper bound on that rectangle is
\[
\begin{aligned}
&k_0(G(a)+C(a))d\,\mathcal E\!\left(\frac{c}{b+\varepsilon_0}\right)\\
&+\frac{1+1/(25a)}{2.999}
 \left(3+\frac{(d+26/25)^4}{4a^2}\right)
 \mathcal E\!\left(\frac{2c}{b+\varepsilon_0}\right). \tag{4}
\end{aligned}
\]
This uses the upper \(x\)-endpoint in the polynomial factor and the lower endpoint in the exponential factor, so it remains valid throughout each rectangle. The interval \(t\le21/20\) covers the entire negative-\(\Pi\) region by (1).

Here is the verifier. All assertions and maxima use `Fraction`; floats are used only to display the maxima.

```python
from fractions import Fraction as F

def E(y):
    return (1 + y/F(128))**(-128)

k0 = F(2601, 14995)
eps = F(63, 1250)

t = F(1, 50)
small = (k0*13/t**2*E(F(1250,88))
         + F(3000,2999)*(3 + F(51,25)**4/(4*t**2))*E(F(1250,44)))
assert small < F(1,100)

mid_max = F(-1)
for i in range(20, 100):
    a, b = F(i,1000), F(i+1,1000)
    r = k0*(8*(1-a)/a**2 + F(19,50)/(a*(a+F(1,50))**2))*E(1/(b+eps))
    u = (1+F(1,25)/a)/F(2999,1000)*(3+F(51,25)**4/(4*a**2))*E(2/(b+eps))
    mid_max = max(mid_max, r+u)
assert mid_max < F(1,3)

top_max, top_arg, count = F(-1), None, 0
for i in range(20, 210):
    a, b = F(i,200), F(i+1,200)
    g = 8*(1-a)/a**2 if a < 1 else F(0)
    corr = F(19,50)/(a*(a+F(1,50))**2)
    pref = (1+F(1,25)/a)/F(2999,1000)
    for j in range(200):
        c, d = 1+F(j,200), 1+F(j+1,200)
        r = k0*(g+corr)*d*E(c/(b+eps))
        u = pref*(3+(d+F(26,25))**4/(4*a**2))*E(2*c/(b+eps))
        count += 1
        if r+u > top_max:
            top_max, top_arg = r+u, (i,j)
assert count == 38000 and top_max < F(9,10)

print(float(small), float(mid_max), float(top_max), top_arg)
```

The exact comparisons give \(0.007963<1/100\), a middle-bin maximum \(0.267909<1/3\), and a top-bin maximum \(0.884427<9/10\), at \(t\in[63/200,8/25]\), \(x\in[1,201/200]\). Hence \(RB+UB^2<1\) throughout the long-suffix region.

The main agent re-ran this verifier.  All assertions pass, and the maxima
are `0.00796`, `0.268` and `0.8844`.

### Appendix to FM44 (three even factors, verbatim from FM-L123i; inserted after FM-CHK11)
The three-factor inequality

Write \(H_b=[b,0]\), \(V=[1,0]\), \(A=[1,1]\), and \(H_0=[0,0]\). Set
\[
Q=H_2-3A+2V+6H_0.
\]
Then the target is \(\Delta(W)=\langle W,Q\rangle\ge0\). This is the \(m=r=2\) reduction recorded in [FM39](</home/yang/q3adjoint/ginibre_q3/SU2_FUNDAMENTAL_MINUS_REDUCTION_2026_09_26.md:2609>).

**Proposition.** For all \(a,b,c\ge0\),
\[
\Delta(H_a\otimes H_b\otimes H_c)\ge0.
\]

**Proof.** Use the supplied multiplicity-free formula for \(H_a\otimes H_b\). Put \(s=a+b\), \(d=|a-b|\), and \(m=\min(a,b)\). Its constituents are
\[
[s-2j-k,k],\qquad 0\le j\le m,\quad 0\le k\le m-j.
\]
At fixed second coordinate \(k\), the first coordinates are exactly
\[
d+k,\ d+k+2,\ \ldots,\ s-k. \tag{1}
\]
Call this support \(R(a,b)\).

By Frobenius reciprocity, \(m_{H_aH_bH_c}(\lambda)\) is the number of constituents in \(R(a,b)\) that occur in \(\lambda\otimes H_c\). The supplied product formula gives \(H_2\otimes H_c\) and \(V\otimes H_c\). For \(A\otimes H_c\), the \(SO(5)\) vector rule gives
\[
A\otimes H_c=
\begin{cases}
A,&c=0,\\
[2,1]\oplus A\oplus V,&c=1,\\
[c+1,1]\oplus[c,1]\oplus[c-1,1]\oplus[c,0],&c\ge2.
\end{cases}
\]
Here the vector rule is
\[
[x,y]\otimes V
=\bigoplus_{\substack{\mu=(x\pm1,y),(x,y\pm1)\\\mu\ \mathrm{dominant}}}[\mu]
\ \oplus\ [x,y]\quad(y>0);
\]
it follows by multiplying the \(B_2\) Weyl alternant by the vector weights \(0,\pm e_1,\pm e_2\). Also \(V\otimes V=H_2\oplus A\oplus H_0\), which gives the displayed \(A\otimes H_c\) formula. Thus all the multiplicities used below are \(0\) or \(1\).

For \(c=0\), \(\Delta(H_aH_b)\) is read directly from \(R(a,b)\):
\[
\begin{array}{c|ccccc}
d & (s,d)=(0,0) & d=0,\ s\ge2 & d=1 & d=2 & d\ge3\\ \hline
\Delta(H_aH_b) & 6&4&2&1&0
\end{array}
\]
This also agrees with the accepted two-factor result [FM-L123f](</home/yang/q3adjoint/ginibre_q3/SU2_FUNDAMENTAL_MINUS_REDUCTION_2026_09_26.md:2754>).

For \(c=1\), the small products combine to
\[
Q\otimes H_1=[3,0]-2[2,1]+2[2,0]+4[1,0]-[1,1]+2[0,0].
\]
Testing these six weights against (1) gives
\[
\begin{array}{c|c}
\text{condition on }(s,d)&\Delta(H_aH_bH_1)\\ \hline
s=0&2\\
s\ge2\text{ even},\ d=0&3\\
s\text{ even},\ d=2&2\\
s\text{ even},\ d\ge4&0\\
s=d=1&4\\
s\ge3\text{ odd},\ d=1&3\\
s\text{ odd},\ d=3&1\\
s\text{ odd},\ d\ge5&0
\end{array}
\]
These cases exhaust the possibilities because \(d\equiv s\pmod 2\).

For \(c\ge2\), the explicit products are
\[
\begin{aligned}
H_2H_c={}&[c+2,0]+[c+1,1]+[c,2]+[c,0]+[c-1,1]+[c-2,0],\\
VH_c={}&[c+1,0]+[c,1]+[c-1,0],\\
AH_c={}&[c+1,1]+[c,1]+[c-1,1]+[c,0].
\end{aligned}
\]
Consequently
\[
\begin{aligned}
\Delta(H_aH_bH_c)={}&
\mathbf1_{[c+2,0]\in R}
+2\mathbf1_{[c+1,0]\in R}
+4\mathbf1_{[c,0]\in R}
+2\mathbf1_{[c-1,0]\in R}\\
&+\mathbf1_{[c-2,0]\in R}
+\mathbf1_{[c,2]\in R}
-2\mathbf1_{[c+1,1]\in R}
-\mathbf1_{[c,1]\in R}
-2\mathbf1_{[c-1,1]\in R}.
\end{aligned} \tag{2}
\]
Let \(\eta_-,\eta_0,\eta_+\) indicate whether \([c-1,1],[c,1],[c+1,1]\), respectively, lie in \(R\). From (1), a level-\(1\) point \([x,1]\in R\) implies \([x+1,0]\in R\); if \(x=c+1\), it also implies \([c,0]\in R\). Thus:

- \(\eta_-\) forces the positive \(4\mathbf1_{[c,0]\in R}\), covering its negative weight \(2\).
- \(\eta_0\) forces the positive \(2\mathbf1_{[c+1,0]\in R}\), covering its negative weight \(1\).
- \(\eta_+\) forces both \([c+2,0]\) and \([c,0]\), giving positive weight \(1+4=5\), covering its negative weight \(2\).

The level-\(1\) first coordinates have fixed parity, so \(\eta_0\) cannot occur together with either \(\eta_-\) or \(\eta_+\). If both \(\eta_-\) and \(\eta_+\) occur, their combined negative weight is \(4\), while the forced positive weight is at least \(5\). This proves (2) is nonnegative; its omitted terms are nonnegative as well.

At the smallest \(c\ge2\) boundary, \(a=b=0,c=2\), the exact margin is \(1\), from \([0,0]\subset H_2H_2\). This completes the proof for all three labels. ∎

### FM-T1-DOM (luna_max_saturn): domination screen for products of symmetric powers

- **Screen.**  All products of 2, 3 and 4 symmetric powers with labels
  `<= 6` and `a <= 8` (1,827 cases): no negative `L(M;a)`.
- **Kills.**
  - `m_M(p,q) <= m_M(p+q,0)` fails.  For `M = (Sym^1)^(x)3`,
    `m(2,1) = 2 > m(3,0) = 1`.
  - Pure-row-only payment fails.  For `M = Sym^2 (x) Sym^4 (x) Sym^6`,
    `a = 8`: the pure rows give 1866, the negative demand is 1879, and the
    positive non-pure terms (2015) are needed.  `L = 2002`.
- **Surviving conjecture (row interval).**
  `m_M(p,q) <= sum_(j=0)^q m_M(p+q-2j, 0)`.  The maximum ratio seen is
  `1294/1389` (four 6s, `(p,q) = (9,1)`).  It is proved for two factors.
- **Surviving scalar companion.**
  `-Lambda(p,q;a) <= sum_(j<=q) Lambda(p+q-2j, 0; a)` whenever
  `Lambda(p,q;a) < 0`: 29 negative cases, `a <= 8`, minimum slack 2.
  Combining the two overcounts shared pure rows, so a non-overcounting
  payment scheme is still needed.


### Corollary FM46 (two symmetric powers, smaller one at most 11; FM-T1-2L10, luna_max_saturn; accepted)

For `t = 10, 11` the same shape holds:

```text
L = binom(n+2v,n+v)^2 (s+1)(s+2)(s+3) F_t Pi_t / D_t,
F_10 = (a+1)(a+2),   F_11 = (a+1)^2 (a+2)(a+3),
Pi_t = sum_(j=0)^22 A_(t,j)(n) v^j.
```

Tables:

```text
t = 10

A_10,0 = (n+1)(n+2)(n+3)^2(n+4)^2(n+5)^2(n+6)^2(n+7)^2(n+8)^2
         (n+9)^2(n+10)^2(n+11)(n+12)(n+13)(n+14)(n^2-3n+6)
A_10,1 = 2P(2232906180599808000,8922873989178777600,16011495601992529920,17365643077047012864,12992376687094698624,7331303850642031104,3379718006445162432,1369059995384504976,508235808714674984,171751897978666936,50982219206210584,12850276294625267,2694004846690344,464834018638038,65625085191056,7542677787493,701103389464,52178845100,3062256408,138457389,4650584,109222,1600,11)
A_10,2 = 2P(5284167471139737600,19321509375727119360,32136021371708455680,32589879584476483008,22887646435109265696,12096347078762079360,5172748027366300176,1914605225414474140,641432614394888954,195001348335210358,52226802028582017,11909650134688750,2257521698463067,350695034638052,44237940892014,4493514778928,363688814560,23100602582,1125335597,40555974,1018123,15888,116)
A_10,3 = 4P(3814240942028021760,12806446617894352896,19756819526566950144,18690102760288756320,12262507649610244512,6037593693918881328,2387853747049717408,809093391412014470,246002677392361342,67637872585567760,16381005803876872,3374607937111319,575620109029052,79882630989598,8905661108948,787781266876,54427593530,2870331028,111484068,3003855,50124,390)
A_10,4 = P(14712549598298810880,46162248527827748736,66713343439164342336,59060685084209774400,36141256696704911520,16510503150483707696,6016859925705088576,1864333394131473410,515292766350300211,128402510452669808,28121906103411130,5217577174111068,795731448083296,97674775944052,9489955097156,717071135374,41158724003,1731616588,50312794,901620,7506)
A_10,5 = 2P(5242541110761778560,15466903579028854656,20911986441117732288,17232632907112929792,9768767163605916032,4114223513935068472,1375316947788052112,389048437750716889,97821378633698516,22109935042613128,4372942141610424,727166432366520,98291462730064,10532124994436,875258512240,55005361685,2521196524,79421156,1536120,13746)
A_10,6 = 2P(2964115790356913280,8131135799466366720,10176354266271032160,7736327681443914088,4034188942680571148,1557796221207966400,475666194966430520,122502593587561390,27978778844444921,5726904608557588,1019004035443708,150738064292084,17841067272076,1639075467488,113502304030,5704845590,196051711,4116332,39806)
A_10,7 = 4P(675232112774060160,1702207266260155776,1957461708274011664,1366690448041320140,653079698701975516,230167866716709254,63854456448974884,14899888864264121,3081154875102430,569433407766053,90639525150740,11802940939840,1202763805106,92363408251,5127670208,193618549,4442828,46720)
A_10,8 = P(983493895285339680,2282983052945828376,2415311937589546076,1547432525173185762,674854323815337875,215585137353309624,53923898390349422,11335771360760568,2116278695874578,351952532017500,49712431046788,5612279529336,480514747702,29695002684,1242152872,31397286,361711)
A_10,9 = 6P(48764575172510120,104228260407149912,101471519185340164,59374047315139801,23432414903233236,6719042895569452,1502811699634460,283304347902460,47634868660000,7089499599484,877039285524,83884963640,5815033128,271999632,7644600,97347)
A_10,10 = 2P(36245826988488480,71669258222984730,63641883225366195,33621924963899958,11871512918477895,3025716763617716,600702409387946,101069383788288,15234500698656,2008752827040,213275572056,16684834708,881833486,27867300,396546)
A_10,11 = 4P(3816496538020710,6880559393747196,5508198790685964,2599334475713625,813882615879842,182766034523179,31898759117976,4751181076584,637697069034,73467216831,6491660272,391933976,14113554,227487)
A_10,12 = 2P(1401829211333340,2209630477665678,1573333100729973,658603523288044,181322285087484,35358279599748,5333992963212,697801034514,83197795437,8203556028,571041022,23798528,442052)
A_10,13 = 52P(7255607087310,11146747555206,7121395420306,2634828933377,626809056192,103342973628,13236078894,1533062208,162181116,13024038,638594,13981)
A_10,14 = 52P(1103674238130,1248456660690,704403254135,224298247788,44536817652,6017881314,658127541,69090198,6234081,363990,9691)
A_10,15 = 3432P(1075302070,1783819312,867796388,230761720,36766170,4060722,390728,36295,2500,85)
A_10,16 = 429P(1751829520,1126083504,460041184,95700962,12303531,1111824,90786,6766,323)
A_10,17 = 858P(6050640,36518064,10967140,1949391,192116,12572,696,61)
A_10,18 = 286P(18464760,4428750,1408185,174206,11543,-56,52)
A_10,19 = 572P(-220830,119052,19268,2185,-110,5)
A_10,20 = 286P(55380,5106,1571,-78,1)
A_10,21 = 3432P(-65,31,-1)
A_10,22 = 17160

t = 11

A_11,0 = (n+2)(n+3)(n+4)^2(n+5)^2(n+6)^2(n+7)^2(n+8)^2(n+9)^2
         (n+10)^2(n+11)^2(n+12)(n+13)(n+14)(n+15)(n^2-3n+6)
A_11,1 = P(136307694995619840000,425877966398158848000,607726409777543577600,530463083232669419520,324220765908803466240,153355960055128415040,61618107436751620032,22581997114392077136,7665477391148096752,2329098972543148396,607302424786766228,132301854049770647,23785174896818829,3511195896594117,424411176851059,41856567674938,3348078108626,215103522518,10932427770,429407739,12567553,257929,3311,20)
A_11,2 = P(191396615742839808000,582388737926331340800,810434500141453079040,689749057881993432960,408874876439881648320,184755460563863741760,69040260643530833952,22873554085261636024,6935451223913962660,1893918179319133804,448688997902675868,89385634091884675,14701552003343387,1976761737849557,215854092012981,19007641571918,1336509470186,73929074894,3144247966,99180423,2184647,29985,193)
A_11,3 = 2P(97465074733719244800,276129807150765219840,359533648287608018688,287724361035868624512,160612887343758338640,68053124651094092664,23592064330615381220,7155454747804361894,1970572206533943950,488852252106440093,105511155516280736,19155230074401395,2860380835962949,346550301594829,33717527096518,2605508569006,157512684634,7282617097,248504118,5895033,86787,597)
A_11,4 = P(138968462220643745280,369499842697177855872,453712547857584567936,342551845636823770784,179832063306736910080,71176518393300586152,22817316701485603028,6331367643469049002,1585038284148136753,357103090190596448,70005692435906142,11511882851130334,1546599934637082,166786055345750,14229564096266,944830377548,47720588183,1769577050,45399460,719700,5310)
A_11,5 = P(74389227571160928000,186168069401594300544,215239498700255479648,152620067021421279920,74895239242289598056,27533626640612574388,8137346325269779542,2067290935794451709,471832272634018637,96700617559250895,17185763183891851,2544448282146571,304455844640831,28803193101479,2112199000301,117200596076,4748603120,132422834,2270958,18048)
A_11,6 = P(31592320586915495040,74873235775945853728,81144925981886052656,53562693832697604472,24339029077415107060,8243407576547206430,2233214770679012717,517992084592310783,107687210368334828,20045718085795805,3215206036251882,424936918379389,44687969280786,3638922760823,223232046899,9944943728,303245558,5656778,48670)
A_11,7 = 2P(5586413508649501920,12231125702990740992,12273460838757614584,7489528115912453408,3135720085833361852,974030871601350487,240849690193556563,50837066714232304,9604949141676709,1618865951445339,232818830770502,27153114307058,2465430236269,168324393622,8304732493,278829220,5694580,53362)
A_11,8 = P(3262930197564828000,6610589142624191144,6102236293201037408,3419388088631933710,1309637478631896414,370047669499086346,82839038461215363,15817028547953842,2706678611481210,411317281053784,52564519617079,5323727579886,406989859159,22415464380,835432324,18824892,193507)
A_11,9 = P(760843987745510160,1437441703765379880,1233600482090160966,636387284651483549,222040068484769845,56627899804588801,11389993990244725,1958929633449373,302959482485749,41274342643547,4624038866211,397187773671,24621183639,1028262759,25800225,293460)
A_11,10 = P(162211088631244920,273895250064878406,212323335466695021,98840942138832299,30918321109982040,7026702105562051,1258768014579316,194145545837155,27027179953456,3266613945377,314058321532,22067291925,1043762148,29487159,375099)
A_11,11 = 2P(12162806315175690,20697851777811594,14946388635035486,6294939703748003,1755855365121096,352906126416268,55878196881057,7684158117151,956375061946,100918319732,8059665411,436751571,14079858,202995)
A_11,12 = P(4616417505464820,5862754085305546,3621169955975861,1338119702704954,327418546910367,57297825758472,7897532414671,960978125198,106302842441,9563308508,600168169,22434750,372879)
A_11,13 = 26P(13087594188270,22733771078226,13415940373552,4454837938707,949486448313,141585148068,16713550840,1809128218,177014464,12922355,570237,11172)
A_11,14 = 26P(2940719759250,2395962845254,1143380521465,321735997721,57267550972,7022356138,711338572,70029730,5837944,308905,7355)
A_11,15 = 52P(15129886410,81026764026,37976491764,9209590552,1325070002,132937181,11990207,1057118,66953,2029)
A_11,16 = 26P(29666584740,13132673238,4484556525,825366211,94449328,7813735,621673,43948,1846)
A_11,17 = 13P(-1179266040,1188160464,352110502,56672433,5018917,322023,18445,1344)
A_11,18 = 13P(301808160,53871574,13496861,1437271,91316,343,371)
A_11,19 = 182P(-589770,133878,20158,2275,-100,5)
A_11,20 = 91P(85260,5750,1595,-74,1)
A_11,21 = 1092P(-105,29,-1)
A_11,22 = 5460
```

**Checks (main agent).**

- *Identities:* both are exact symbolic identities in `(s,v)` (342 s and
  443 s; the coefficient tables were parsed directly from the packet).
- *Coefficient signs (Sturm counts in `n`):* every `A_(t,j)` is positive on
  `[0, inf)` except the following.
  - `A_(10,19)`: root at `n ~ 1.458`, so positive for `n >= 2`.
  - `A_(11,17)`: root at `n ~ 0.785`, so positive for `n >= 1`.
  - `A_(11,19)`: root at `n ~ 2.841`, so positive for `n >= 3`.
  - `A_(t,21)`: the middle coefficient of the final quadratic
    `A_(t,20) + A_(t,21) v + A_(t,22) v^2`.  Its discriminant has no root
    in `[0, inf)` (value `-1037402308800` for `t = 10`, `-156302218800` for
    `t = 11`), so the quadratic is positive.
- *Small `n`:* the whole polynomial `Pi_t(n,v)` in `v` has no root in
  `[0, inf)` and is positive at `v = 0`, for `n = 0, 1` (`t = 10`) and
  `n = 0, 1, 2` (`t = 11`), by Sturm.

**Corollary FM46.**  Q3 holds for minus labels `{s+1, t+1, 1, 1}` with any
number of plus 1s whenever `min(s,t) <= 11`.

The packet also records that the uniform coefficient-shape conjecture fails
at `t = 11` (`A_(11,17)(0) < 0`); the per-`t` certificates still close.
The per-`t` line is paused.  The structural route is FM-T1-DOM.

### Theorem FM47 (the same-sign `{1,3}` family is complete; FM-T2b5, luna_max_venus; accepted)

*Local notation* (added after FM-AUDIT1): in FM47 and its appendix,
`q = m/N` and `beta = b/N` are local ratios, and `dmu = dsigma(x) dsigma(y)`
with `sigma` the semicircle law.  `M(m,r)` is the FM18 normaliser, and
`Phi` is the rectangle ratio `Phi_rect`.  These are unrelated to the later
`q`, `beta_a`, `M(u,v;a)` and `Phi(d,t,a)`.

**Claim.**  For every odd `b >= 7` and integers `0 <= m < r`,
`I(b,2m,r) > 0`.

**Proof** (packet FM-T2b5).  Put `N = m+r`, `q = m/N` and `beta = b/N`.

1. *Negative part.*  On `{H_3 < 0}`, `3U + V + 4(-H_3) = 8` with
   `U = u^2`, `V = v^2`.  Weighted AM--GM gives
   `U^m V^r |H_3|^b <= B = (8m/(3(N+b)))^m (8r/(N+b))^r (2b/(N+b))^b`.
   So the negative part is at most `B`.
2. *Positive rectangle.*  Take
   `R_s = [2-eps, 2-eps/2] x [4s-2-eps/4, 4s-2+eps/4]` in the semicircle
   coordinates `(x, y)`, with `eps = 1/200`; at `(2, 4s-2)`, `u = 4s` and
   `v = 4(1-s)` (clarified after FM-CHK13),
   which lies inside `{H_3 > 0}`.  Its measure is `> c_0 = 3/1228800000`.
   With `s` chosen by regime (`s = q`, `1/100`, `1/2`, `4/5`), the minimum
   integrand over `R_s` is `>= (7/5)^N B`.
3. *Closing.*
   - For `N >= 60`: `c_0 (7/5)^60 > 1`.
   - For `b >= 2N + 24`: `c_0 (5/4)(5/2)^24 > 1`.
   - In both cases `I > 0`.
4. *Finite box.*  `N <= 59`, `b <= 2N + 23`: exact screen of 43,655 cases,
   all positive.  The minimum normalised value is `38247989/53954384`, at
   `(m,r,b) = (6,21,7)`.

**Checks (main agent).**

- Re-derived the AM--GM bound, the rectangle measure
  (`3 eps^(5/2)/2048` with `sqrt(eps) > 1/15`, `pi^2 < 16`), the formula
  for `Phi`, and the first row of the table by hand.
- Checked the four-row per-`N` bound numerically: 200,000 random
  `(q, beta)`, minimum log-margin above `log(7/5)` equal to `0.067`.
- Re-ran the packet's exact Fraction verifier: every assertion passes,
  including the case count, the minimum value and the edge values.

**Theorem FM47.**  Combining FM18, FM20, FM22 and this result: Q3 holds for
every word with labels in `{1,3}` whose 3s all have the same sign.  With
FM28 this extends to any number of 2s on the matching half-quadrant.

### FM-CHK12 (luna_max_mercury): FM42 appendix and FM44 forcing

- **FM42 long-suffix verifier: repaired.**
  - As first printed, the `small` term passed `E(F(88,1250))`.  The
    exponent at `t = 1/50` is `1/(t+eps) = 1250/88`.  With the wrong
    argument the value is 5254.26 and the first assertion fails.  The
    main agent's earlier re-run printed this comparison as `False` and it
    was missed.
  - The appendix now reads `E(F(1250,88))`.  Re-run by the main agent:
    - `small = 0.0079622 < 1/100`;
    - `mid_max = 0.2679085 < 1/3`;
    - `top_max = 0.8844265 < 9/10`, over 38,000 rectangles.
    All three pass.
  - Mercury also checked the supporting interval arguments and accepts
    them after this correction:
    - monotonicity for `t <= 1/10`;
    - the rectangle endpoint choices;
    - `e^(-y) <= (1+y/128)^(-128)`;
    - bin coverage up to `t <= 21/20`.
  - Theorem FM42 stands.
- **FM44 `c >= 2` forcing argument: ACCEPTED.**  Mercury checked:
  - the `B_2` alternant identities for `H_2 H_c`, `V H_c` and `A H_c`;
  - formula (2), with level-zero weights `1,2,4,2,1`;
  - the parity exclusion, and the covering weights `4`, `2` and `5`;
  - the boundary values for `c = 0, 1, 2`.

### Corollary FM48 (three plus labels, support range; FM-T2-PQR3, luna_max_mars; accepted)

Word `[1^-, 1^-, p^+, q^+, s^+, 1^+ x a]` at `r = 1`, with `p >= q >= s`.

**Fusion identity.**  For every permutation `(u,v,w)` of `(p,q,s)`:

```text
T(p,q,s;a) = sum_(j in C(u,v)) T(j,w;a) + sum_(j in C(u,w)) M(j,v;a) + sum_(j in C(v,w)) M(j,u;a)
```

Here `C(u,v) = {|u-v|, |u-v|+2, ..., u+v}`.  The first sum is nonnegative:
use FM42 for `j >= 1` and FM24 for `j = 0`.  Take `(u,v) = (p,q)`, fusing
the two largest labels.  Every correction then has index sum at least
`p+q-s`, and `M(u,v;a) = 0` when `u+v > a+2`.

**Corollary FM48.**  `T(p,q,s;a) >= 0` whenever `p+q-s > a+2`.  Also
`T = 0` when `p+q+s+a` is odd.

**Residual.**  The open case is `p >= q >= s`, `p+q-s <= a+2`, with
`p+q+s+a` even.  It is exactly inequality (3):

```text
X + sum_(j in J1) M(j,q;a) + sum_(j in J2) M(j,p;a) >= 0,   X = sum_(j in C(p,q)) T(j,s;a)
```

Here `J1 = {j in C(p,s) : j+q <= a+2}`, `J2 = {j in C(q,s) : j+p <= a+2}`
and `R = sum_(J1) M(j,q;a) + sum_(J2) M(j,p;a)`.  Negative corrections
already occur on the support boundary, e.g. `(3,2,2;3)`: `X = 11`,
`R = -3`, `T = 8`.

**Screen and checks.**

- *Exact screen:* `2 <= s <= q <= p <= 20`, `3 <= a <= 30`, 8,988 residual
  tuples, all `T >= 0`.
  - `R < 0` in 4,103 cases.
  - Minimum `T = 7`, at `(5,2,2;3)`.
  - Minimum of `T/(X + sum|negative corrections|)` is `7/15`, at
    `(3,3,2;6)`: `X = 183`, `R = -64`.
- *Main agent:* re-ran the packet verifier with identical output.  Checked
  identity (1) independently against a direct evaluation of
  `E[(x-y)^2 S_p S_q S_s (x+y)^a]` over the semicircle product (sympy
  expansion with Catalan moments).  In all 121 cases (`p, q, s <= 5`,
  `a <= 6`) the direct value equals `2T`.

### Theorem FM49 (one arbitrary label of each sign at `r = 1`; FM-T2-PM, luna_max_jupiter, plus main-agent certificate)

**Claim.**  Q3 holds for `[q^-, 1^-, p^+, 1^+ x a]` for all `p, q >= 1`
and `a >= 0`.

**Reduction** (jupiter).
- `q >= p` is FM39, and `p = 1` is FM11.
- For `p > q`, FM5 applied to the weight `(p-1,q)` gives
  `S_p H_q = sum_(j in C(p,q)) H_j - Res V_(p-1,q)`.  So the word value
  is `J = 2 Delta_(p,q)(a)` (FM6/FM7), with
  `Delta_(p,q)(a) = sum_(j in C(p,q)) m_a(j-1,0) - m_a(p-1,q)`.
- Here `m_a` is the multiplicity in `(C^4)^(x)a`.  By Brauer's formula
  (1) and `beta_a(t) = (2t/(a+1)) binom(a+1, (a+1)/2 + t)`:
  - `m_a(A,B) = alpha(u) beta(v) - beta(u) alpha(v)`, with
    `u = (A+B+3)/2` and `v = (A-B+1)/2`;
  - `alpha(t) = beta(t-1) + beta(t) + beta(t+1)`.
- The one-row part telescopes:
  `m_a(A,0) = g(v) - g(v+1)`, where `g(v) = beta(v)^2 - beta(v+1) beta(v-1)`.

**Closed form.**  The word value is 0 unless `a+p+q` is odd (parity); assume
that.  Put `h = (a-p-q+1)/2 >= 0` (otherwise the subtracted
term vanishes), `S = p+q`, `D = p-q`, `N = a+1` and
`rho = binom(N,p+h)/binom(N,h) = prod_(i=1)^q (1 + p/(h+i))`.  Then

```text
Delta = binom(N,h)^2 Psi(rho) / ( N (S+h+1)^2 (S+h+2) (p+h+1) (q+h+1) ),
Psi(x) = (D^2+S+2h)(S+h+1)^2(S+h+2) x^2 - D(p+1)(q+1)(S+2)(S+h+1)(S+2h+1) x
         - h(p+h+1)(q+h+1)(S^2+5S+2h+4).
```

`Psi` has positive leading coefficient and `Psi(0) <= 0`.  So
`Delta >= 0` iff `rho >= x_+`, the positive root.

**Proof of `Psi(rho) >= 0`.**
- *`q = 1, 2`* (jupiter).  `Psi(rho)` equals a positive factor times
  `Q_q(p,h)`, where
  `Q_1 = 12h^2 + (-4p^2 + 4p + 36)h + p^4 + 2p^3 - 5p^2 + 2p + 24` and
  `Q_2 = 12h^2 + (-4p^2 + 4p + 60)h + p^4 + 2p^3 - 9p^2 + 6p + 72`.  Its discriminant in `h` is
  `-16(2p^4 + 8p^3 + 2p^2 - 12p - 9) < 0`.
- *`q >= 3`* (main agent).
  1. *Jensen.*  `i -> log(1 + p/(h+i))` is convex, so
     `rho >= (1+y)^q >= B_3 = sum_(k<=3) binom(q,k) y^k`, with
     `y = 2p/(2h+q+1)`.
  2. Since `B_3 > 0`, it suffices that `Psi(B_3) >= 0`.
  3. Put `q = 3+Q`, `p = 4+Q+P`, `h = H`.  Then
     `36(2h+q+1)^6 Psi(B_3)` is a polynomial in `P, Q, H` with 689
     terms, of total degree 17.  Exactly 14 coefficients are negative,
     at `P^2 Q^j H^7` (`j <= 3`), `P^3 Q^j H^6` (`j <= 5`) and
     `P^3 Q^j H^7` (`j <= 3`).
  4. Each negative term sits at the exponent midpoint of two positive
     terms in the same `Q^j` slice:
     - `(2,7)` between `(0,8)` and `(4,6)`;
     - `(3,6)` between `(1,7)` and `(5,5)`;
     - `(3,7)` between `(1,8)` and `(5,6)`.
     In every case `c_neg^2 <= 4 c_1 c_2`, and no positive term is used
     twice.  So AM--GM gives `Psi(B_3) >= 0` on the whole orthant
     `P, Q, H >= 0`.

**Checks** (main agent; verifier `character_ring_iter/verify_fm49_mixed_sign.py`, runs in 15 s).

- (A) The closed form is an exact symbolic identity in `(p, q, h, rho)`,
  starting from the telescoped `beta` expression.
- (A') The regrouped determinant form, the `alpha`-`beta` relation and
  identity (2) agree with Brauer's formula (1) for every `a <= 24`.
- (B) The `q = 1, 2` factorisations and discriminants hold symbolically.
- (C) The polynomial expansion and the AM--GM pairing, in exact integers.
- Brauer's formula (1) agrees with direct dominant-weight peeling of
  `(C^4)^(x)a` for `a <= 6`.
- `J = 2 Delta` agrees with direct semicircle evaluation of
  `E[D_q D_1 S_p (x+y)^a]` for `p <= 6`, `q < p`, `a <= 7`.
- Jupiter's exact screen: 143,370 triples (`p <= 60`, `a <= 80`), no
  negative value.

**Independent check:** ACCEPTED in FM-CHK14 (mercury).

*FM-T2-PQR4 (mars), addendum to FM48.*
- **Screen maximum.**  The largest ratio is `(-R)/X = 11275/32174 < 0.351`,
  at `(4,4,2;12)`.  Main-agent scan (`p, q, s <= 30`, `a <= 40`):
  - for `s = 2`, the ratio stays near `0.33-0.35` along `p = q`,
    `a ~ p(p-1)`;
  - it decreases with `s` (`s >= 10`: `<= 0.06`).
- **`s = 2` reduction.**  `S_2 = (S_1^2 + D_1^2)/2 - 2` gives
  `T(p,q,2;a) = T(p,q;a+2)/2 + B_2(p,q;a)/2 - 2T(p,q;a)`, where `B_2` is
  the `r = 2` word with two plus labels.
- No uniform bound yet.  Reassigned to the FM49 template.

### Theorem FM50 (`r = 2` core, four even factors; FM-L123k, luna_max_neptune; accepted)

**Claim.**  `Delta(H_a H_b H_c H_d) >= 0` for all `a, b, c, d >= 0`, where
`Delta(M) = <M, Q>` and `Q = [2,0] - 3[1,1] + 2[1,0] + 6[0,0]` (SO(5)).

**Proof** (packet).
1. *Kernel.*  `chi_Q = (X-Y)^2`, with `X = z_1 + 1/z_1` and
   `Y = z_2 + 1/z_2`.  In coordinates `U = x+y`, `V = x-y` its weight
   stencil is `kappa(u) kappa(v)` at `(2u, 2v)`, where `kappa(0) = 2`,
   `kappa(+-1) = -1`, and `kappa = 0` otherwise.
2. *Pair decomposition* (FM44 appendix):
   `H_a H_b = sum_(0 <= j <= i <= m) [(U,V) = (alpha+2i, alpha+2j)]`, with
   `alpha = |a-b|` and `m = min(a,b)`; for the second pair `(c,d)`,
   `beta = |c-d|` and `n = min(c,d)`.  Here `H_a = [a,0]` is the SO(5)
   harmonic of the `Sp(4)'` picture and `Q` is the SO(5) virtual module.
   These are distinct from the original `H_q` and the Sp(4) weights
   `Q_r`.
3. *Brauer--Klimyk matching.*  Match the two pair supports against the
   stencil.  Only the identity and one wall reflection contribute (equal
   parity), or only the coordinate swap (opposite parity).
4. *Sum.*  With `b_t = e_t - e_(t+1)`, the double sums become interval
   overlap counts.

This gives an exact formula.  In the equal-parity case, with
`h = (alpha-beta)/2`:

```text
Delta = ov([h,h+m+1],[0,n+1]) + ov([h,h+m],[0,n]) + 1{h=0, m=n} + 1{alpha=beta=0}(1 + 1{m=n})
```

The opposite-parity case, with `h = (alpha-beta-1)/2`, is:

```text
Delta = ov([h+1,h+m+1],[0,n+1]) + ov([h,h+m+1],[0,n])
```

Here `ov([r,s],[u,v]) = max(0, min(s,v) - max(r,u) + 1)` is the number of
integers in the intersection of two closed intervals (definition added
after FM-CHK14).  Every term is nonnegative.

**Checks (main agent).**
- Verified `chi_Q = (X-Y)^2` from the three SO(5) characters, and the
  stencil conversion.
- The formula equals a direct computation of `Delta` for all
  `a, b, c, d <= 6` (2,401 cases).  The direct computation takes the
  harmonic characters `Sym^j - Sym^(j-2)`, multiplies by `(X-Y)^2` and the
  Weyl denominator, and reads off the constant term.

**Open.**  Five or more factors.  The packet gives `Delta(P H_e)` in terms of
the low multiplicities `p_(x,y)` of a four-factor product `P`, for
`e = 0`, `1` and `>= 2`.

### Theorem FM51 (row-interval domination for all products; FM-T1-DOM2, luna_max_saturn; accepted)

**Claim.**  For every partition `lambda` and every `p >= q >= 0`:

```text
b_lambda(p,q) <= sum_(j=0)^q b_lambda(p+q-2j, 0)
```

Here `b_lambda(p,q)` is the multiplicity of `V_(p,q)` in `S_lambda(C^4)`
restricted to `Sp(4)`.  Hence (*) holds for every polynomial GL(4)-module,
in particular for every product of symmetric powers.

**Proof.**
1. `SL(4) = Spin(6)` contains `Sp(4) = Spin(5)`, the stabiliser of the
   nonisotropic vector `omega` in `Lambda^2 C^4`.
2. `Spin(6) -> Spin(5)` branching is multiplicity-free with interlacing
   `x_1 >= y_1 >= x_2 >= y_2 >= |x_3|` (Goodman--Wallach Thm 8.1.4).  In
   these coordinates:
   - `x = (b + (a+c)/2, (a+c)/2, (c-a)/2)`, with `(a,b,c)` the
     differences of the parts of `lambda`;
   - `y = ((p+q)/2, (p-q)/2)`.
3. Every constituent `V_(p,q)` comes with the pure row `(2x_2, 0)`.  Its
   weight `(x_2, x_2)` interlaces.  Branching also needs the lattice
   condition `y - x_2` integral, i.e. `|lambda| = p+q (mod 2)`.  So
   `j = (p+q-2x_2)/2` is an integer, `0 <= j <= q` by interlacing, and
   `2x_2 = p+q-2j`.

**Main-agent check.**  The interlacing branching reproduces
`dim S_lambda(C^4)` exactly for all `lambda_1 <= 11`.

**Flow formulation** (packet).
- Demands `m_M(lambda)(-Lambda(lambda;a))` are routed to capacities
  `m_M(mu) Lambda(mu;a)` along row-interval overlaps.  `L >= 0` follows
  from the capacitated Hall condition (2).
- Exact screen: 1,827 product-suffix cases and 26,067 demand subsets, with
  no violation.  The tightest normalised slack is `2/65`.  There is an
  explicit flow for `(2,4,6)`, `a = 8`.
- (*) alone cannot give Hall.  The vector `V_(1,1) + V_(2,0)` satisfies (*)
  but has `L = -2`; it is not a product of symmetric powers.

### Two-row reduction of the two-symmetric-powers family (main agent)

**Branching.**  For a two-row `lambda = (d+t, t)`,
`S_lambda(C^4)|_(Sp(4)) = sum_(q=0)^t V_(q+d, q)`, a diagonal string (from
FM51's interlacing).  Pieri gives

```text
Sym^s W (x) Sym^t W = sum_(j=0)^t S_(s+t-j, j) W,     s >= t
```

so

```text
L(s,t,a) = sum_(j=0)^t Phi(s+t-2j, j, a),   Phi(d,t,a) = sum_(q=0)^t Lambda(q+d, q; a).
```

**Conjecture TR (two-row positivity)** [now Theorem FM54; the statements below are historical].  `Phi(d,t,a) >= 0` for all
`d, t, a`.  It implies the whole two-symmetric-powers family, i.e. all
`s, t`, including the open `t >= 12` range.
- Exact screen: 53,900 cases (`a <= 48`, all `d, t`), no negative value.
- The first differences are `Phi(d,0,a) = Lambda(d,0;a) >= 0` (FM19).

**Telescoping.**  Consider the 24 Weyl-binomial terms of
`Lambda(q+d, q; a)`.
- The factor `binom(a, (a-X+Y)/2)` does not depend on `q`.
- The Weyl pairs `(u,v)` and `(-v,-u)` carry opposite signs and opposite
  `sigma = u+v`.

Hence

```text
Phi(d,t,a) = G(d,0,a) - G(d,t+1,a),   G(d,s,a) = sum_(q>=s) Lambda(q+d,q;a)
           = sum_(pairs) eta_l det(w) binom(a,(a-d-1+delta)/2) * sum_(i=0)^(sigma-1) binom(a,(a+d+3-sigma)/2+s+i)
```

This is a window sum of at most 5 binomials, with `delta = u-v` and
`sigma in {1,3,5}`.  So two-row positivity is `G(d,s,a) <= G(d,0,a)` for
all `s >= 1`.  In the screen, `max G(d,s,a)/G(d,0,a) = 0.946`, at `s = 1`,
where the difference is exactly the pure row.

**Three-row modules fail individually.**  For example
`L(S_(2,1,1) W; 0) = -2`; there are 953 negative three-row cases with
`lambda_1 <= 15`, `a <= 15`.  So three or more factors need the Kostka
weights of the product, as FM51's obstruction also indicates.

**Three factors via Pieri (main agent, exact screen).**  With
`k_1 >= k_2 >= k_3`:

```text
h_k1 h_k2 h_k3 = sum_(j=0)^(k_2) h_k3 s_(k1+k2-j, j)
```

Write `nu = (d+t, t)`.
- `L(h_k s_nu; a)` has 145 negatives in 9,477 cases when `k` is
  unrestricted.  The first is `L(h_2 s_(1,1); 0) = -1`.
- Restricted to the range that occurs, `k <= t + floor(d/2)`: 28,067 cases
  and 24 negatives.  **All** of them are `nu = (t,t)` with `k = t`, i.e.
  the all-equal product `h_t^3` at `j = t`.
- The full three-factor sums `L(h_k1 h_k2 h_k3; a)`: 2,475 cases
  (`k_i <= 8`, `a <= 14`), no negative value.

So, apart from one term, fusing the two largest factors gives termwise
positivity for three factors.  The exception is `h_t s_(t,t)` inside
`h_t^3`, which must be paid for by its siblings `j < t`.

### Appendix to FM47 (verbatim from FM-T2b5, luna_max_venus; inserted after FM-CHK13)

Coordinates: `x, y` are the two semicircle variables, `u = x+y`, `v = x-y`, `U = u^2`, `V = v^2`, `H = x^2+xy+y^2-2` (`= H_3`), `I(b,2m,r) = int U^m V^r H^b dmu`.

#### A.2 Bound the negative part and construct a positive rectangle

On H < 0, set Z = −H = 2 − (3U+V)/4, so 3U + V + 4Z = 8. Weighted AM–GM gives the pointwise bound

UᵐVʳZᵇ ≤ B,

where

B = (8m/(3(N+b)))ᵐ (8r/(N+b))ʳ (2b/(N+b))ᵇ.

For m=0 the first factor is interpreted as 1. Since μ is a probability measure, the integral of the negative part has magnitude at most B.

Set ε=1/200. For 0<s≤4/5 let yₛ=4s−2 and define

Rₛ = [2−ε, 2−ε/2] × [yₛ−ε/4, yₛ+ε/4].

At (2,yₛ), the values of u=x+y, v=x−y and H are 4s, 4(1−s), and hₛ=16s²−8s+2. On Rₛ,

u ≥ 4s−1/160, v ≥ 4(1−s)−1/160, H ≥ hₛ−1/25.

Every rectangle used below lies in [−2,2]². The semicircle density is √(4−x²)/(2π). On the x interval, √(4−x²)>√ε; on all chosen y intervals, |y|≤1569/800 and 4−y²>9/64. Using π<4 and √ε>1/15 gives

μ(Rₛ) > 3ε⁵ᐟ²/2048 > c₀, c₀ := 3/1,228,800,000.

Let Pₛ=(4s)²ᵐ(4(1−s))²ʳhₛᵇ. Direct simplification gives

Pₛ/B = Φ(q,β;s)ᴺ,

Φ(q,β;s) = 2(1+β)3ᑫ (s²/q)ᑫ ((1−s)²/(1−q))¹⁻ᑫ (hₛ(1+β)/(2β))ᵝ,

with the q=0 term understood by continuity.

The following bounds compare the rectangle’s minimum integrand to Pₛ. The final column is a lower bound for the ratio to B per N.

| Parameter range | s | Boundary ratio Φ | Rectangle factor | Combined lower bound |
|---|---:|---:|---:|---:|
| β≤1, q≥1/50 | q | 3/2 | (99/100)²·24/25 | 88209/62500 > 7/5 |
| β≤1, q≤1/50 | 1/100 | 441/250 | >19/20 | 8379/5000 > 7/5 |
| 1≤β≤2 | 1/2 | 2 | (319/320)²·(49/50)² | 244328161/128000000 > 7/5 |
| β≥2 | 4/5 | (6/25)(73/25)ᵝ | (99/100)²(99/100)ᵝ | ≥1440747/781250 > 7/5 |

Here are the checks behind the table.

- For s=q and β≤1, hq=16(q−1/4)²+1≥1 and 3ᑫqᑫ(1−q)¹⁻ᑫ≥3/4; the latter has its minimum at q=1/4. [Corrected after FM-CHK14: the raw relative losses need not be ≤1/100 — at q=1/50 the raw u-loss is 5/64 — but the powered bounds hold: (u_min/(4q))^(2q) ≥ 99/100 (log loss ≤ 1/295) and (v_min/(4(1−q)))^(2(1−q)) ≥ 99/100.] Its relative loss in H is at most 1/25.
- For q≤1/50 and s=1/100, hₛ=1201/625. Concavity of q↦q log((3/10000)/q) shows (3s²/q)ᑫ≥(3/200)¹ᐟ⁵⁰>9/10; the other q factor is at least 9801/10000. The rectangle factors are bounded below by (134/135)(599/600)²(1176/1201)>19/20.
- At s=1/2, ∂q log Φ = log(3s²(1−q)/((1−s)²q))>0 for 0<q≤1/2. At q=0, Φ=(1+β)(1+1/β)ᵝ/2, which increases for β≥1 and equals 2 at β=1. The local factors follow from u₀=v₀=2 and hₛ=2.
- At s=4/5 the same derivative is positive. At q=0, Φ=(2/25)(1+β)((73/25)(1+1/β))ᵝ, giving the stated bound for β≥2. On the rectangle, each of u, v and H is at least 99/100 of its boundary value.

Consequently, for every parameter pair, the rectangle’s minimum integrand divided by B is at least (7/5)ᴺ. Therefore, whenever N≥60,

∫Rₛ u²ᵐv²ʳHᵇ dμ / B ≥ c₀(7/5)ᴺ ≥ c₀(7/5)⁶⁰ > 1.

The rectangle lies in H>0, so its contribution is positive. Subtracting the negative-part bound B proves I(b,2m,r)>0 for all N≥60.

For β≥2 the last row gives the stronger bound

(minimum integrand on R₄⁄₅)/B ≥ (1/5)ᴺ(5/2)ᵇ.

Thus if b≥2N+24, then

c₀(1/5)ᴺ(5/2)ᵇ ≥ c₀(5/4)ᴺ(5/2)²⁴ ≥ c₀(5/4)(5/2)²⁴ > 1.

This proves positivity for that range at every N.

#### A.3 Exact finite screen

The only remaining parameters satisfy

1≤N≤59, 7≤b≤2N+23, b odd, 0≤m<r=N−m.

The Fraction verifier checks all 43,655 cases. It finds no nonpositive value; the minimum normalized value is 38247989/53954384 at (m,r,b)=(6,21,7). Since M(6,21)=138816056788291168, the corresponding integral is 98405998168043878.

Small boundary values, computed from the same exact formula, are

- I(7,0,1)=1222;
- I(7,0,4)=1820 and I(7,2,3)=3052;
- I(7,0,7)=3586258 and I(7,2,6)=18538.

The verifier uses only integer and Fraction arithmetic:

```python
from fractions import Fraction as F
from math import factorial

def L_values(m, r, B):
    N = m + r
    Am, Ar = [1], [1]
    for i in range(B):
        Am.append(Am[-1] * (2*m+2*i+1) * (2*m+2*i+3))
        Ar.append(Ar[-1] * (2*r+2*i+1) * (2*r+2*i+3))

    den = [1]
    for n in range(B):
        den.append(den[-1] * (N+2+n) * (N+3+n))

    c = []
    for n in range(B+1):
        c.append(sum((
            F(3**i * Am[i] * Ar[n-i],
              factorial(i) * factorial(n-i) * den[n])
            for i in range(n+1)
        ), F(0)))

    return [
        factorial(b) * sum((
            c[n] * F((-2)**(b-n), factorial(b-n))
            for n in range(b+1)
        ), F(0))
        for b in range(B+1)
    ]

def M(m, r):
    return F(
        2 * factorial(2*m) * factorial(2*m+1)
          * factorial(2*r) * factorial(2*r+1),
        factorial(m)**2 * factorial(r)**2
          * factorial(m+r+1) * factorial(m+r+2)
    )

c0 = F(3, 1228800000)
assert c0 * F(7, 5)**60 > 1
assert c0 * F(5, 4) * F(5, 2)**24 > 1
assert F(59, 64) > F(99, 100)**25
assert F(3, 200) > F(9, 10)**50
assert F(4) - F(1569, 800)**2 > F(9, 64)
assert F(3, 2) * F(99, 100)**2 * F(24, 25) > F(7, 5)
assert F(134, 135) * F(599, 600)**2 * F(1176, 1201) > F(19, 20)
assert F(441, 250) * F(19, 20) > F(7, 5)
assert 2 * F(319, 320)**2 * F(49, 50)**2 > F(7, 5)
assert F(6, 25) * F(99, 100)**2 * F(14, 5)**2 > F(7, 5)

count = 0
best_value, best_arg = None, None
for N in range(1, 60):
    B = 2*N + 23
    for m in range((N+1)//2):
        r = N - m
        values = L_values(m, r, B)
        for b in range(7, B+1, 2):
            value = values[b]
            assert value > 0
            count += 1
            if best_value is None or value < best_value:
                best_value, best_arg = value, (m, r, b)
    print(f"N={N}/59, checked={count}", flush=True)

assert count == 43655
assert best_value == F(38247989, 53954384)
assert best_arg == (6, 21, 7)
assert M(6, 21) * best_value == 98405998168043878

edges = {
    (0, 1, 7): 1222,
    (0, 4, 7): 1820,
    (1, 3, 7): 3052,
    (0, 7, 7): 3586258,
    (1, 6, 7): 18538,
}
for (m, r, b), expected in edges.items():
    assert M(m, r) * L_values(m, r, b)[b] == expected
```


### FM-CHK13 (luna_max_mercury): FM47 and FM46

- **FM46: ACCEPTED.**  Mercury re-derived the exact Sturm isolations of the
  exceptional coefficients:
  - `A_(10,19)` root in `(4556/3125, 589/404)`;
  - `A_(11,17)` root in `(2053/2615, 2232/2843)`;
  - `A_(11,19)` root in `(6041/2126, 1452/511)`.
  Also checked: the final-quadratic discriminants `D_10(n)` and `D_11(n)`,
  negative on `n >= 0`; the small-`n` full-polynomial checks; and the
  domain `v >= 0`.
- **FM47: reported DEFECT; resolved as a reading issue.**
  - Mercury read `R_s` in `(u,v)` coordinates.  Then `s = 1/2` would put
    `v = 0` in the rectangle, where the integrand vanishes.
  - The packet's rectangle is in `(x,y)` coordinates, and there `u` and `v`
    stay near `4s` and `4(1-s)`.  The note's statement has been clarified.
  - The missing details are now in the note: the regime table, `Phi`, the
    rectangle factors and the verifier, as the appendix above.
  - Mercury confirmed the AM--GM bound, the measure bound, the closing
    constants and the 43,655-case count.
  - The appendix was re-checked in FM-CHK14 (bounds valid, one sentence corrected).

*Correction to the three-factor Pieri observation (main agent, larger
range).*  The rule "only `nu = (t,t)`, `k = t` is negative" holds only in
the first screen (`t <= 12`, `a <= 16`).  For `h_t^3` with `t <= 14` and
`a <= 30`:
- the sibling term `h_t s_(t+1,t-1)` is also negative in several cases
  from `t = 11` on (e.g. `t = 11`, `a = 15`);
- pairing the square term with its nearest sibling fails in 26 of 434
  cases.

Still, all 434 totals `L(h_t^3; a)` are positive, and the negative Pieri
terms are at most `0.35%` of the positive ones (worst `1/290`, at
`t = 6`, `a = 6`).  So the fused-largest-pair grouping is a
small-perturbation structure, not termwise positivity.

### Packets FM-T2-PQR5 (jupiter), FM-T1-2L11 (venus), and two main-agent kills

- **FM-T2-PQR5 (jupiter), three plus labels.**
  - Two exact formulas:
    - `2T = sum_(I subset {p,q,s}) sum_r binom(a,r) [B_(r+2)(I) B_(a-r)(I^c) - 2 B_(r+1)(I) B_(a-r+1)(I^c) + B_r(I) B_(a-r+2)(I^c)]`,
      with `B_d(I)` built from SU(2) fusion and ballot numbers;
    - `T = sum_lambda c_lambda m_a(lambda)`, with
      `hat S_t = V_(t,0) - V_(t-1,1) + V_(t-2,0)`.
  - No uniform bound.  The residual of FM48 remains open.
- **FM-T1-2L11 (venus), two symmetric powers at `t >= 12`.**
  - Proved the boundary `a = s-t` (`v = 0`):
    `L(t+n, t, n) = (n^2-3n+6)/2 > 0` for all `t, n >= 1`, via explicit
    chamber-walk counts of the three surviving constituents.
  - `M_2 >= 0` fails along this boundary (`M_2 = -2`), and at the first
    open point (`t = 12`, `n = 1`, `a = 3`: `M_2 = -3`, `L = 6`).
  - Open at the time: `t >= 12`, `n >= 1`, `v >= 1` [closed by FM54].
- **Kill K22 (main): per-`mu` positivity at `r = 1`.**
  - By Schur--Weyl, `T = sum_mu f^mu P(mu)`, with
    `P(mu) = sum_(lambda interlacing mu) c_lambda`.
  - For a single `hat S_p`, `P(mu) >= 0` holds: the FM51 pure-row
    companions give `x_2 in {p/2, (p-2)/2}`.
  - It fails for products: `hat S_2 hat S_2` at `mu = (2,1,1)` gives
    `P = -1`.  11 of 21 two-label products and 24 of 35 three-label
    products have negative `P(mu)`.
- **Kill K23 (main): inward-flow certificate uniform in `a`.**
  - `T(a) = sum kappa(gamma) b(k_1) b(k_2)`, with
    `b(k) = binom(a, a/2+k)` symmetric and log-concave.  A flow of positive
    `kappa`-mass inward in weak majorisation would give `T >= 0` for all
    `a` at once.
  - It fails already for the word `[1^-,1^-,1^+,1^+ x a]`
    (flow 3 < demand 4).  Log-concavity alone is too weak.
- **Observation (main).**  The `r = 2` kernel
  `Lambda(lambda;a) = <V_lambda (x) W^(x)a, V_(2,0) - 3V_(1,1) + 5V_(0,0)>`
  has character `(X-Y)^2`, the same as FM50's `Q`.  By Jacobi--Trudi,
  `Phi(d,t,a) = L(d+t,t,a) - L(d+t+1,t-1,a)`.  So Conjecture TR says the
  two-power value decreases as the labels spread apart with their sum
  fixed.

- **FM-T2-QQP (mars), `[q^-, q'^-, p^+, 1^+ x a]`.**
  - Reduction: the only open cases are `p > q >= q' >= 2`, `a >= 5`,
    `p <= q+q'+a`, with even parity.  The rest is covered by FM49, FM11,
    FM39 and Cor. 23A9ZZ10.
  - Exact formula for `J` as Catalan moments.  Screen: 2,197 cases
    (`p <= 12`, `a <= 24`), all positive; the minimum is 2, at
    `(9,2,2;5)`.
  - No closed form yet.  Restarted in a fresh thread.

**Corollary FM52 (`q' = 2` slice; FM-T2-QQP2, luna_max_mars; accepted).**
FM1 gives `D_2 = D_1 S_1`, i.e. `chi_2(x) - chi_2(y) = (x-y)(x+y)`.  So
`[q^-, 2^-, p^+, 1^+ x a] = [q^-, 1^-, p^+, 1^+ x (a+1)]`, and FM49 gives
Q3 for all `p, q, a`.
- The packet's exact check compares direct Catalan moments with
  `2 Delta_(p,q)(a+1)`: 685 cases, all agree.
- The two-minus family `[q^-, q'^-, p^+, 1^+ x a]` stays open for
  `p > q >= q' >= 3`, `a >= 5`.
- The next slice does not reduce the same way.  `chi_3 = x^3 - 2x` gives
  `H_3 = x^2 + xy + y^2 - 2 = (S_1^2 + S_2 - S_0)/2`, with `S_0 = 2`.
  (Corrected: a first version omitted the `-S_0` term.)  So `q' = 3` equals
  `(1/2)[q^-,1^-,p^+,1^+ x (a+2)] + (1/2)[q^-,1^-,p^+,2^+,1^+ x a] - [q^-,1^-,p^+,1^+ x a]`.
  This mixes signs, so the direct `S_p H_q H_(q')` route is preferable.

### TR closed form and the `r = 2` kernel as a Laplacian of the `r = 1` kernel (main agent, with FM-T1-TR2 venus and FM-T2-R2PM jupiter)

Let `M(u,v;a) = E[(x-y)^2 (x+y)^a U_u(x) U_v(y)]`.  This is the FM30
kernel, with closed form
`M = beta_(e-1) beta_(d-1) 16 (N-2e+2)(u+1)(v+1) P_M / (...)`.  Let
`tau_d(a)` be the FM24 one-label value.

1. **Integral identity** (venus, checked).
   - `Phi(d,t,a) = (1/2) E[(x-y)^2 (x+y)^a (D_(d+t+1) D_(t+1) - D_(d+t+2) D_t)]`.
   - It follows from the bialternant
     `(x-y)^2 s_(d+t,t)(z,1/z,w,1/w) = D_(d+t+1) D_(t+1) - D_(d+t+2) D_t`
     together with `U_m U_n - U_(m+1) U_(n-1) = U_(m-n)`.
2. **Three-term formula.**
   `Phi(d,t,a) = tau_d(a) - M(d+t+1, t+1; a) + M(d+t+2, t; a)`.
   Checked exactly in 756 cases.  `tau_d = (1/2) E[(x-y)^2 (x+y)^a S_d]`
   exactly.  The FM30 closed form also holds at `q = 1` (108 cases) and
   `p = q` (78 cases).
3. **Laplacian** (convention `U_(-1) = 0`, i.e. `M(u,-1;a) = 0`).  Differencing in `t` gives
   `Lambda(P,Q;a) = M(P,Q) + M(P+2,Q) - M(P+1,Q+1) - M(P+1,Q-1)`
   (checked in 2,310 cases).  So every `r = 2` H-only word value is a
   combination of `r = 1` kernels.
4. **TR as one ratio condition.**  Put `D = d/2`, `e = D+t+1`, `N = a+2`
   and `beta_x = binom(N, N/2+x)`.  The two `M` terms share `e`, so
   `Phi / beta_D^2 = R_A - rho K`, where:
   - `rho = beta_(e-1)/beta_D = prod_(i=1)^t (H+1-i)/(N-H+i)`, with
     `H = N/2 - D`;
   - `R_A = tau_d/beta_D^2` and `K` are explicit rational functions.

   Hence TR holds iff `rho K <= R_A`.
5. **Jensen.**
   - `log((N-H+i)/(H+1-i))` is convex in `i`, so
     `1/rho >= (1+y)^t`, with `y = (4D+2t)/(2H+1-t)`.  This is exact at
     `t = 1`.
   - With the truncation `B_k = sum_(j<=k) binom(t,j) y^j` and
     `N = 2D+2t+2+2V`, `t = 1+T`, the numerator of `R_A B_k - K` is a
     polynomial in `(D,T,V)`.  It has 192/286/424 terms, of which
     42/56/82 are negative, for `k = 1, 2, 3`.  So no immediate
     certificate.
6. **Asymptotically tight.**
   - The worst ratio `rho K / R_A` occurs at `d = 0`, `t = 1`, and tends
     to 1 as `V -> infinity`: `1 - ratio ~ 178/V^2` (0.868 at `V = 30`,
     0.99982 at `V = 1000`).  Similar limits hold for other `d, t`.
   - So TR is tight at leading order as `a -> infinity`.  A certificate
     must match the leading asymptotics exactly: an expansion in `1/V`,
     then a region split.

*Remark (main agent): TR at `t = 1` as a certificate check of the method.*

At `t = 1`, `rho = H/(N-H+1)` exactly, and

```text
Phi(d,1,a)/beta_D^2 = 4(D+2)(2D+1)(2D+3) Q(D,V) / [(V+3)(D+V+2)(2D+V+3)^2(2D+V+4)^2(2D+V+5)(2D+V+6)(2D+2V+3)],
Q = 8D^5 + 12D^4V + 40D^4 + 4D^3V + 14D^3 - 36D^2V^2 - 127D^2V - 114D^2
    + 24DV^2 + 91DV + 84D + 45V^3 + 285V^2 + 570V + 360.
```

`Q > 0` on the orthant, by AM--GM on exponent midpoints:
- `D^2V^2` from `(4/5)12 D^4V` and `(4/5)45 V^3`: `4 * 48/5 * 36 >= 36^2`;
- `D^2V` from `(3/4)40 D^4` and `(1/2)285 V^2`: `17100 >= 127^2`;
- `D^2` from `(1/4)40 D^4` and `360`: `14400 >= 114^2`.

This re-proves the `t = 1` two-row pieces, which are already covered by
FM26.  The open content of TR is uniformity in `t`, assigned as FM-T1-TR3.
The rational identity uses the FM30 closed form at `q = 1` and, for
`d = 0`, at `p = q`.  These edge cases were checked exactly (108 and 78
cases) but are not separately proved symbolically.

### FM-CHK14 (luna_max_mercury): FM49, FM51, FM47 appendix, FM50

- **FM49: ACCEPTED.**  Mercury checked:
  - the reduction and the Brauer determinant (spot checks for `a <= 9`);
  - the telescoping and the closed form;
  - the `q = 1, 2` factorisations;
  - Jensen, and an independent expansion of the `q >= 3` polynomial
    (689 terms, the same 14 negatives, each AM--GM pair valid, no reuse).
  - Note: the verifier prints some Booleans rather than asserting them.
- **FM51: ACCEPTED.**  Stabiliser, `D_3` labels, interlacing, lattice
  parity, and the small cases.
- **FM50: ACCEPTED.**  Independent `B_2` constant-term check for all labels
  `<= 4`.  The `ov` convention is now stated in the note.
- **FM47 appendix: bounds valid; one sentence corrected.**  The first-row
  justification speaks of raw relative losses `<= 1/100` in `u` and `v`.
  At `q = 1/50` the raw `u` loss is `5/64`.  What the proof needs, and what
  holds, are the powered bounds:
  - `(u_min/(4q))^(2q) >= 99/100`, with log loss at most `1/295`;
  - `(v_min/(4(1-q)))^(2(1-q)) >= 99/100`.

  The table factors are unaffected.  Mercury replayed the finite box
  exactly (43,655 cases), with the same minimum and edge values.

### Theorem FM53 (two arbitrary minus labels, one arbitrary plus label; FM-T2-QQP3, luna_max_mars; accepted)

**Claim.**  Q3 holds for `[q^-, q'^-, p^+, 1^+ x a]` for all
`p, q, q' >= 1` and `a >= 0`.

**Proof** (mars).
1. Let `E = Sym^(q-1) W (x) Sym^(q'-1) W (x) W^(x)a`, a polynomial
   GL(4)-module.
2. By FM6/FM7 and self-duality,
   `J/2 = <hat S_p, E>_(Sp(4)) = b_E(p,0) - b_E(p-1,1) + b_E(p-2,0)`,
   since `hat S_p = V_(p,0) - V_(p-1,1) + V_(p-2,0)`.
3. FM51 at `(p-1,1)` gives `b_E(p-1,1) <= b_E(p,0) + b_E(p-2,0)`.
   This is for `p >= 2`.  For `p = 1` (added after FM-CHK15),
   `hat S_1 = V_(1,0)` and `J/2 = b_E(1,0) >= 0` directly.

(For label `0`: `hat S_0 = 2 V_(0,0)`, so `<hat S_0, E> = 2 dim E^(Sp(4)) >= 0`.  This
covers the `j = 0` terms used later.)  Per GL(4) constituent the slack
`B_p(lambda)` is 0 or 1.  This subsumes
FM49 and FM52, and the `r = 1` single-plus-label cases of FM39.

**Checks (main agent).**
- Direct semicircle integral `= 2 sum_lambda c_lambda(E) B_p(lambda)`, and
  every `B_p(lambda) >= 0`: 1,701 cases (`q <= 6`, `q' <= q`, `p <= 9`,
  `a <= 8`, both parities).
- The interlacing branching must include the lattice condition
  `|lambda| = A+B (mod 2)`.  Without it, spurious negative slacks appear.
  The packet's verifier omits the condition but only runs parity-matched
  cases, where it is automatic.
- The corrected branching reproduces `dim S_lambda` for all
  `lambda_1 <= 11`.

### Unified formulation and Conjecture MP_2 (main agent)

**Unified formulation.**
- For symmetric `f`, the Sp(4) Weyl density relative to
  `sigma (x) sigma` is `(x-y)^2/2`.  So a word with `2r` minus labels
  `q_i`, plus labels `p_j >= 2` and `a` plus 1s has value
  `2 <hat S_(p_1) ... hat S_(p_k) (x) E, Q_r>_(Sp(4))`, where:
  - `E = (x)_i Sym^(q_i - 1) W (x) W^(x)a` is a polynomial GL(4)-module;
  - `chi_(Q_r) = (x-y)^(2(r-1))`, with `Q_1 = 1` and
    `Q_2 = V_(2,0) - 3V_(1,1) + 5V_(0,0)`.
- `chi_(Q_2) = (2-X')(2-Y')`, where `X' = zw + 1/(zw)` and
  `Y' = z/w + w/z` are the SO(5) torus coordinates.  This is the product
  over the short roots of `|1 - e^alpha|^2`.
- Status:
  - `r = 1`: `k = 0` is FM11 and `k = 1` is FM53; `k = 2` with `E = W^a`
    is FM42.  The general `E` is open for `k >= 2`.
  - `r = 2`, `k = 0`: next item.

**The `r = 2`, `k = 0` case as a Kostka inequality.**
- Interlacing gives `w(lambda) = <S_lambda W, Q_2>` in `{5, 2, 1, -2, 0}`.
  With `(a, b, c)` the successive part differences of `lambda`:
  - `w = 5` for `a = b = c = 0`;
  - `w = 2` for `a = c = 0`, `b >= 1`;
  - `w = 1` for `(a,c) in {(2,0), (0,2)}`;
  - `w = -2` for `a = c = 1`;
  - `w = 0` otherwise.

  Checked against `L(S_lambda;0)` in 364 cases.
- So Q3 for every word with four arbitrary minus labels and plus 1s is
  `sum_lambda K_(lambda,kappa) w(lambda) >= 0`, with
  `kappa = (q_1-1, ..., q_4-1, 1^a)`.

**Correction (main agent).**  This is not a new conjecture.  With
`w(lambda) = m_(2,0) - 3m_(1,1) + 5m_(0,0)` evaluated on `S_lambda`, the
statement below is exactly the `r = 2` H-only part of Conjecture FM3, via
Lemma FM12 (same weight vector).  Its multigraph form is the fifth-pass
target `3 m_(1,1) <= m_(2,0) + 5 m_(0,0)`, in the 3-nesting-free model.
The label "MP_2" is kept only as shorthand.

**Conjecture MP_2 (= FM3, `r = 2`, H-only).**  `F_2 = sum_(l(lambda)<=4) w(lambda) s_lambda` is
monomial-positive, i.e. `<h_kappa, F_2> >= 0` for every composition
`kappa`.
- Exact screen: all 915 partitions of size `<= 16`, no negative value
  (464 zeros).
- It contains the whole two-symmetric-powers family,
  `kappa = (s, t, 1^a)`, and every `{q_1, ..., q_4}` minus-label family
  with plus 1s.
- A combinatorial proof would inject SSYT of the `w = -2` shapes into those
  of the positive shapes.

**Conjecture MP_r (all `r`).**
`F_r = sum_(l(lambda)<=4) <S_lambda W, Q_r> s_lambda` is monomial-positive.
It is equivalent to Q3 for every word whose plus labels are all 1 (any
`2r` minus labels).
- `Q_3 = V_(4,0) - 5V_(3,1) + 10V_(2,2) + 14V_(2,0) - 35V_(1,1) + 35V_(0,0)`,
  from `chi = (x-y)^4`.
- Exact screens: `r = 2`, all partitions of size `<= 16` (915);
  `r = 3`, all partitions of size `<= 14` (508).  No negative value.
- Integral form: `<h_kappa, F_r> = (1/2) E_(sigma (x) sigma)[(X-Y)^(2r) prod_i sum_(a+b=kappa_i) U_a(X) U_b(Y)]`,
  using `det(1 - x g)^(-1) = (sum_k U_k(X) x^k)(sum_l U_l(Y) x^l)` on Sp(4).

**Conjecture MP = Conjecture FM3, restated (main agent; see the
correction above).**  The weight vectors of `Q_2` and `Q_3` coincide with
Lemma FM12's `r = 2, 3` vectors.  The following is FM3 in GL(4)/Sp(4)
language, with new screens.  For every `r >= 1` and every
multiset of plus labels `p_1, ..., p_k >= 2`, the symmetric function

```text
F = sum_(l(lambda)<=4) <S_lambda W, hat S_(p_1) ... hat S_(p_k) (x) Q_r>_(Sp(4)) s_lambda
```

is monomial-positive: `<h_kappa, F> >= 0` for every composition `kappa`.
- Words have `kappa = (q_1-1, ..., q_(2r)-1, 1^a)`.  So MP implies Q3 for
  all words.
- For general `kappa` (more than `2r` parts `> 1`), MP is the H-only
  statement `E[(x-y)^(2r) prod_(i=1)^n H_(q_i) prod_j S_(p_j)] >= 0` with
  `n` arbitrary.  The note already uses such products (FM44, FM50).
- It is not Conjecture HQ (K20): HQ was a `(T,S)` cone in the `Sp(4)'`
  picture.
- **Proved cases.**  `k = 0`, `r = 1` (trivially); `k = 1`, `r = 1` for
  every `kappa`, by FM51 as in FM53 (per `lambda`).
- **Exact screens, all partitions `kappa` of size `<= 11`**, no negative
  value:
  - `r = 2` with plus `(2)`, `(3)` or `(4)`;
  - `r = 1` with plus `(2,2)`, `(3,2)`, `(3,3)`, `(4,2)` or `(2,2,2)`;
  - `r = 3` with plus `(2)`.

  With `k = 0`: `r = 2` up to size 16, and `r = 3` up to size 14.
- **Per-`lambda` positivity fails** once `k >= 2` (`r = 1`, K22) or
  `r >= 2` (e.g. `w(2,1,1) = -2`).  So a proof must use the `h_kappa`
  product structure: Kostka numbers, Pieri, or RSK.

- **Additional packets.**
  - FM-T1-TR3 (venus): TR proved at `t = 1` for all `d`, independently of
    the main-agent remark.  It includes the edge case `M(2,2)` at `d = 0`
    via the new identity `M(p,1;a) = tau_p(a+1) - tau_(p+1)(a) - tau_(p-1)(a)`,
    which also makes the `q = 1` edge of FM30 rigorous through FM24.  TR
    was then open for `t >= 2` [closed by FM54].
  - FM-T2-R2PM2 (jupiter): normalised the row-interval inequality to
    `A_delta - rho^2 B_E >= rho |Psi|`, in one binomial ratio `rho`, with
    a Jensen bound for `delta >= 1/2`.  No certificate: `rho -> 1` as
    `a -> infinity`, so fixed-margin bounds cannot close it.

- **Correction.**  "MP" and "MP_r" are FM3 restated in Kostka form; see
  Lemma FM12 and the fifth pass.  They are not new conjectures.  The new
  results today toward FM3 are:
  - FM51: row-interval domination for every polynomial GL(4)-module.  It
    makes the fifth-pass `m_(1,1) <= m_(2,0)`-type bounds unconditional,
    in the weaker form `m_(1,1) <= m_(2,0) + m_(0,0)`.
  - FM53: FM3 for `r = 1` with exactly one `hat S_p` factor and arbitrary
    `H_q` factors.  Since FM51 and the product decomposition work for every
    `kappa`, this holds for all such products, not only words.

### Theorem FM54 (Conjecture TR; the two-symmetric-powers family is complete; FM-T1-TR, luna_max_saturn; accepted)

**Claim.**  `Phi(d,t,a) >= 0` for all `d, t, a >= 0`.  Hence Q3 holds for
every word with minus labels `{s+1, t+1, 1, 1}` and any number of plus 1s,
for all `s, t`.  This subsumes FM26, FM35, FM38, FM41, FM43, FM46 and the
`v = 0` boundary result, and closes the open range `t >= 12`.

**Proof** (saturn).  Put `b_j = binom(a,j)`, `n = (a-d)/2` (same
parity), `alpha = A_0`, `beta = B_0 - 3A_0`, `gamma = A_0 - B_0 + C_0`.
1. *Window formula.*
   `G(d,s,a) = A sum_(i<5) b_(d+n-1+s+i) + B sum_(i<3) b_(d+n+s+i) + C b_(d+n+1+s)`,
   with `A, B, C` explicit combinations of `b_(n-3..n+2)`.
2. *Basis change.*  Rewriting
   `1+z+...+z^4 = (1+z)^4 - 3z(1+z)^2 + z^2` and using `b_(a-j) = b_j`
   gives `G(d,s,a) = b_n F(n+1-s)`, where
   `F(x) = alpha binom(a+4,x) + beta binom(a+2,x-1) + gamma binom(a,x-2)`.
3. *Closed forms.*
   - `alpha = (d+1) delta / ((n+1)(d+n+1)(d+n+2))`, with
     `delta = d^2 - d - 6n - 6`.
   - `G(d,0,a) = b_n^2 (d+1)(d+2n+1)(d+2n+2) R_0 / ((n+1)^2(n+2)(d+n+1)^2(d+n+2)^2(d+n+3)) > 0`.
     `R_0` has discriminant `-16(2d^4 + 8d^3 + 2d^2 - 12d - 9)` in `n`,
     negative for `d >= 2`.  For `d = 0, 1`, `R_0 = 12(n+1)(n+2)`
     (corrected after FM-CHK-TR).
4. *Endpoint.*  `F(n+1) >= F(0) = alpha`.  Use
   `R_0 - delta^2 = 2(d+2n+1)(delta + d^2 + 2d)`, the box inequality for
   the three factors, and `b_n >= (n+1)(n+2)/2` when `n >= 2`; for
   `n = 1`, `R_0 - 6 delta > 0`.
5. *Shape.*  For `1 <= x <= n`, the sign of `F(x+1) - F(x)` is the sign of
   an explicit quadratic `Psi` in `y = n+1-x`.  `Psi` is positive at
   `y = 1`, and increasing (if `epsilon > 0`) or concave (if
   `epsilon <= 0`).  So `F` has no interior maximum on `[1, n+1]`.  With
   `F(1) - F(0) = (d+1)H/(...)` and
   `Psi(x=1) = (d+2n+3)(H + 3(8 - delta))`, this gives
   `max_(0<=x<=n+1) F = F(n+1)`.
6. Hence `Phi = G(d,0,a) - G(d,t+1,a) >= 0`.  The support and parity cases
   are separate.

**Checks (main agent).**
- Every displayed polynomial identity verified symbolically (sympy):
  - `alpha`, the difference formula and the `P`--`Psi` factorisation;
  - `Psi(n)` and its discriminant; `G(d,0,a)`; `R_0 - delta^2`;
  - `F(1) - F(0)`, `Psi(x=1)`, the `H + 3(8-delta)` identity and the
    linear-coefficient identity; the `n = 1` identity.
- The inequality steps re-derived by hand, including the `F(1)` versus
  `F(0)` case split and the monotonicity or concavity case analysis.
- The window identity `G = b_n F(n+1-s)` and the maximum property
  (6): exact, 13,299 cases (`a <= 60`).
- The packet verifier (`character_ring_iter/verify_fm54_two_row.py`):
  `PASS`, 45,521 triples (`a <= 48`).
- Independent check: ACCEPTED in FM-CHK-TR (venus), after two repairs; the full proof is in the Appendix to FM54.

### Three-term formula for every Schur module (main agent)

For `lambda = (l_1,l_2,l_3,0)` with part differences `(A,B,C)`, the `r = 2`
value with `a` plus 1s is

```text
L(S_lambda W; a) = M(A, C) - M(A+B+1, B+C+1) + M(A+B+C+2, B),
```

where `M(u,v;a) = E[(x-y)^2 (x+y)^a U_u(x) U_v(y)]` is the FM30 kernel,
and `M(u,0) = tau_u`.
- *Derivation.*  Laplace-expand the `4 x 4` bialternant of
  `s_lambda(z,1/z,w,1/w)` along the rows `(z,1/z) | (w,1/w)`.  Each
  `2 x 2` minor is `(z-1/z) U_(l_j-l_k-1)(X)`, and the cross Vandermonde
  is `(X-Y)^2`.
- Checked exactly in 2,145 cases (`lambda_1 <= 8`, `a <= 12`).  TR is the
  two-row case.
- Consequently the `r = 2` H-only part of FM3 is
  the target inequality (OPEN)
  `sum_lambda K_(lambda,kappa) [M(A,C) - M(A+B+1,B+C+1) + M(A+B+C+2,B)] >= 0`
  (Kostka weights).  Here `M` is the defining integral.  Closed forms are:
  - FM30 for `u > v >= 2` in the interior;
  - `M(u,0) = tau_u` (FM24) for `v = 0`;
  - `M(p,1;a) = tau_p(a+1) - tau_(p+1)(a) - tau_(p-1)(a)` (FM-T1-TR3) for
    `v = 1`;
  - equal labels `u = v`: the FM30 formula at `d = 0`, checked in 78
    cases but not separately proved.

- **Kill K24 (FM-MP-A, venus): Littlewood-factor positivity fails.**
  - At `r = 1`, `k = 2`, `p = q = 2`: `f_00 = 3` and `f_11 = 2`.  The
    quotient by `(1 - t_1 t_2)^(-1)` then has `[t_1 t_2] = -1`.  The
    target coefficients themselves are positive.
  - The stable restriction rule also needs King's modification at rank 4:
    `det W` restricts to the trivial module, whereas the unmodified rule
    gives dimension 6.
  - **Exact reformulation.**  `hat S_q = h_q + 2h_(q-2) - s_(q-1,1)`, so
    the `r = 1`, `k = 2` part of FM3 is
    `A_p((q-1,1);kappa) <= A_p((q);kappa) + 2 A_p((q-2);kappa)`,
    where `A_p(mu;kappa) = <hat S_p, h_kappa s_mu> = sum_lambda c_lambda(h_kappa s_mu) B_p(lambda)`
    and `B_p(lambda) >= 0` by FM51.

## Status ledger (2026-09-27, after FM54)

Every entry was checked by the main agent.  "(chk)" marks an independent
checker pass.

**Proved infinite Q3 families.**

- **Labels in `{1,2}`:** all words.  Repo Cor. 24B17B / Prop. 24B3.
- **Two minus labels (`r = 1`):**
  - **Two arbitrary minus labels, at most one plus label `>= 2`, any
    number of plus 1s: complete** (Theorem FM53, via FM51 (chk)).  This
    subsumes FM24, FM49 (chk) and FM52.
  - It extends to FM3's whole `r = 1`, one-`hat S` sector: any number of
    `H_q` factors.
  - Two plus labels `p, q`, minus labels 1: all `p, q, a` (FM42 (chk)).
  - Three plus labels, minus labels 1: `p+q-s > a+2` (FM48), plus an exact
    screen.
  - `<= 7` factors (Cor. 23A9ZZ10); all plus labels 1 (FM11); absorption
    (FM39); FM21.
- **Four minus labels (`r = 2`):**
  - **`{s+1, t+1, 1, 1}` with plus 1s: complete for all `s, t`**
    (Theorem FM54, TR; independently checked, FM-CHK-TR).
  - `{q, 1, 1, 1}` with plus 1s: FM19, and the twist.
  - Four `Sp(4)'` even factors: FM44 (chk), FM50 (chk: the all-label
    derivation in the Appendix to FM50 was independently ACCEPTED in
    FM-CHK17), and the two-factor
    cases FM45.
- **Labels `{1,3}` and `{1,2,3}`:** same-sign 3s complete (FM47 (chk));
  FM18, FM20, FM22, FM25, FM28, FM29, FM32, FM34.
- **One arbitrary label at every `r`:** labels `<= 5` (FM36, FM37); label
  `>= (N-2)/3` (FM27); `min(a,t) <= 3` (FM33).  Open: the residual branches
  (paused).

**Structural results.**
- FM51: row-interval domination for every polynomial GL(4)-module
  (Spin(6) -> Spin(5)) (chk).
- The unified form of FM3: `<hat S_(p_1) ... hat S_(p_k) (x) E, Q_r>_(Sp(4)) >= 0`,
  with `E` a product of symmetric powers and `chi_(Q_r) = (x-y)^(2(r-1))`.
  FM3 is equivalent to monomial positivity of Kostka-weighted sums
  ("MP", which is FM3 restated).
- Three-term formula
  `L(S_lambda;a) = M(A,C) - M(A+B+1,B+C+1) + M(A+B+C+2,B)` for every
  Schur module.  The `r = 2` kernel `Lambda` is a discrete Laplacian of
  the `r = 1` kernel `M`.

**Open core.**
- FM3 at `r = 2`, H-only, three or more factors (saturn: three factors;
  mars: general, via multigraphs).
- `r = 1` with two or more `hat S` factors and arbitrary `H`'s (the
  reformulation after K24).
- `r = 2` with one `hat S` (jupiter).
- Five `Sp(4)'` even factors (neptune).
- `r >= 3` in general.

**Kills this session:** K20--K24 (HQ cone, per-`mu` positivity at
`r = 1`, inward flow, Littlewood factor).

### FM-CHK15 (luna_max_mercury): FM53, the unified formulation, the TR formulas

- **FM53: ACCEPTED after repair.**
  - The pairing normalisation checks.
  - `hat S_p = V_(p,0) - V_(p-1,1) + V_(p-2,0)` verified independently via
    `chi_(r,s) = (U_(r+1)(X) U_s(Y) - U_s(X) U_(r+1)(Y))/(X-Y)`.
  - `B_p(lambda)` is 0 or 1.
  - The `p = 1` case is now stated separately.
- **Unified formulation and MP_2 reformulation: ACCEPTED.**  Checked the
  Weyl density `(x-y)^2/2`, `chi_(Q_2) = (X-Y)^2 = (2-X')(2-Y')`, and the
  `w(lambda)` table with lattice parity.  MP_2 itself remains conjectural.
- **TR formulas and the `t = 1` certificate: ACCEPTED.**
  - The three-term and Laplacian identities hold, with `U_(-1) = 0` stated.
  - Independent symbolic route for `t = 1`:
    `Phi(d,1,a) = L(d+1,1,a) - L(d+2,0,a)`, via FM26 and FM19.
  - AM--GM margins: `432/5`, `971`, `1404`.
  - `V < 0` is handled by degree support.

- **FM-L123l (neptune), five `Sp(4)'` even factors: exact formula, no proof.**
  - *SO(5) -> SO(4) branching.*  `Delta(prod_i H_(b_i)) = sum_(k <= b) D(k)`
    (sum over `0 <= k_i <= b_i` for each `i`),
    where `c_j` is the multiplicity of `E_j` in `E_(k_1) (x) ... (x) E_(k_n)`
    (SU(2)), and
    - `D = (3c_0 - c_2)(2c_0 + c_2 - c_4)` if `sum k_i` is even;
    - `D = -(2c_1 - c_3)^2` if `sum k_i` is odd.
  - The derivation is valid for every `n`.
  - Checked against the SO(5) Weyl alternant: 252 five-label multisets.
    Screen: 792 multisets through label 7, no negative value.
  - Integral form:
    `Delta = (1/2) int int prod_i K_(b_i)(x,y) (4-x^2)(4-y^2)(x-y)^2 dmu dmu`,
    with `K_b(x,y) = sum_(k<=b) E_k(x) E_k(y)`, a Christoffel--Darboux
    kernel.
  - *Obstruction.*  One-label fibre sums can be negative (labels
    `(1,0,0,0,2)`: fibre sum `-2`, total `2`).  So a fibrewise telescope
    cannot prove it.

- **FM-MP-B (jupiter), `r = 2` with one `hat S_p`: exact weights, no proof.**
  - Explicit branching:
    `S_lambda|_(Sp(4)) = sum_(i<=m, j<=beta) V_(M+i+j, m+j-i)`, with
    `m = min(alpha,gamma)`, `M = max(alpha,gamma)`.
  - The stencil `hat S_p Q_2` in closed form (10 terms for `p >= 4`).
  - A table of the Schur weights `w_p(lambda)`.  Negative (`= -2`) exactly
    on `r in {0,2}` with `beta = 1`, or `r in {4,6}` with `m = 1`, where
    `r = alpha + gamma - p + 4`.
  - The family is `sum_lambda f^(lambda/(s)) w_p(lambda)` (skew SYT
    counts).  Screen: `p <= 8`, `s+a <= 30`, no negative value.
  - Also: `hat S_p = h_p + 2h_(p-2) - s_(p-1,1)` on Sp(4), so this family
    is a three-factor H-only term plus a two-row correction.

### Corollary FM55 (three symmetric factors: `a = 0`, a unit factor, support; FM-T1-3F, luna_max_saturn; accepted)

For `M = Sym^s W (x) Sym^t W (x) Sym^u W (x) W^(x)a`, `r = 2`
(minus labels `{s+1, t+1, u+1, 1}` with plus 1s), `L >= 0` holds in the
following cases.

- **Parity:** `L = 0` if `s+t+u+a` is odd.
- **Support:** with `K` the largest part of `(s,t,u,1^a)` and `R` the sum of
  the other parts, `L = 0` if
  `K - R >= 3`, and `L = m_20 >= 0` if `K - R = 2`.  This uses the FM14
  weight bound.
- **Unit factor:** `min(s,t,u) <= 1` reduces to FM54.
- **The whole `a = 0` slice.**  With at most five vectors, the Sp(4)
  invariant ring of three factors plus the gadgets is freely generated
  (the first Pfaffian relation needs six).  So
  `I(k) = 1_T(k)`, where `T` is the set of even-sum triples satisfying the
  triangle inequalities.  Write `B(k) = 1_T(k)`,
  `D(k) = sum_i B(k - 2e_i)` and `C(k) = sum_(i<j) B(k - e_i - e_j)`
  (indicators are 0 when a coordinate is negative).  The gadget formula
  `4 I(k+(2)) - 3 I(k+(1,1)) + 8 I(k)` becomes `L = 5B + D - 2C`, with
  indicator sums `B`, `D`, `C`.  A short case analysis gives `L >= 0`.

**Checks (main agent).**
- `5B + D - 2C` equals the direct integral for all labels `<= 12`
  (455 triples).
- The packet's screen: 7,735 cases, and 2,238 in the open region, all
  nonnegative; the Pieri/kernel sum matches the direct integral
  everywhere.

**Open** (as first recorded; now `a >= 2` after the `a = 1` slice below): `a >= 1`, `s >= t >= u >= 2`, `s <= t+u+a`, even parity.  There
the fifth-pass 3-nesting constraint is active: `3 + a` factors plus two
gadget vectors make six or more vectors.

- **FM-L123m (neptune): compound-kernel identity and a PSD kill.**
  - *Identity.*  `Delta(H_(b_1)...H_(b_n)) = sum_k (M^k_(b1+1,b2+1) M^k_(b1,b2) - M^k_(b1+1,b2) M^k_(b1,b2+1))`.
    This is a sum of `2 x 2` minors of the moment matrices
    `M^k_(ij) = int E_i E_j prod E_(k_t) (4-x^2) dmu`.
  - *Kill K25.*  Positive semidefiniteness alone cannot give `Delta >= 0`:
    the rank-one kernel `xy` pairs to `-4`.  Single minors are negative
    already for three factors (`H_1^3`: terms `4, -1`).
  - Five-factor screen: 792 multisets, no negative value.  No proof.
- **TR at `r = 3` (main agent, exact).**
  - `Phi_3(d,t,a) = sum_(q<=t) <V_(q+d,q) (x) W^a, Q_3>` is nonnegative in
    1,017 cases (`d <= 14`, `t <= 8`, `a <= 14`), except
    `Phi_3(1,1,1) = -2`.
  - Three-row modules at `r = 3` fail in 392 of 2,420 cases.
  - So the `r = 3` two-symmetric-powers family (minus labels
    `{s+1,t+1,1,1,1,1}`) needs FM54's method plus one small grouping.

### FM-CHK-TR (luna_max_venus, independent checker): Theorem FM54 ACCEPTED

- Window formula and Weyl pairing: accepted.  The coefficients `A, B, C`
  and the window starts were re-derived.
- Basis change: accepted.  The note's `A_0, B_0, C_0` are `A, B, C`
  divided by `b_n`.
- Closed forms and endpoint: accepted after two repairs, both now in the
  note:
  - the discriminant statement for `R_0` (it is `144` at `d = 0, 1`,
    where `R_0 = 12(n+1)(n+2)`);
  - the `n = 0` endpoint: `delta > 0` forces `d >= 4`, and
    `F(1) - F(0) = (d^3 - 2d^2 - 5d + 30)/(2(d+3)) > 0`.  Saturn's
    section 5 treats `n = 0` separately as well.
- Difference factorisation and maximum: accepted.  Venus's `Psi(y)` equals
  saturn's.  An explicit `x = 1` increment
  `F(2) - F(1) = n(d+1)(d+2n+1) Psi(n) / (2(n+1)(n+2)(d+n+2)(d+n+3)(d+2n+3))`.
  Also `Psi(1) = 60(n+1)(n+2)` for `d = 0, 1`.
- Support, parity, and the passage to `L(s,t,a)`: accepted.
- Verifier re-run: 45,521 triples pass.

**FM55, `a = 1` slice (FM-T1-3F2, luna_max_saturn; accepted).**
`L(s,t,u;1) >= 0` for all `s, t, u`.
- With six vector slots the invariant ring is the polynomial ring modulo a
  single `6 x 6` Pfaffian of multidegree `(1^6)`.  So the gadget counts are
  free counts minus one shifted indicator, and
  `L(s,t,u;1) = sum_i A(k-e_i) - 2 sum_i B(k-e_i) + 3 B(k-(1,1,1))`, with
  `A = 5B + D - 2C` (the `a = 0` value).
- `A(k') >= 1` for every triple `k' in T`, and a two-case analysis (`s <= t+u-1`,
  `s = t+u+1`) gives `L >= 0`.
- Main-agent check: the formula equals the direct integral for all
  `s >= t >= u >= 0` with `s <= 10`, `s <= t+u+1` and even parity
  (100 cases).  The packet's screen
  for `a >= 2` has 2,118 cases, all positive (minimum 3).

Open: `a >= 2`, `s >= t >= u >= 2`.

### Round notes (neptune FM-T1-TR3r, venus FM-T2-PQS, mars FM-T1-MP2)

- **`r = 3` two symmetric powers (neptune).**
  - `Phi_3(d,t,a) = M_2(d,0) - M_2(d+t+1,t+1) + M_2(d+t+2,t)`, with
    `M_2(u,v;a) = E[(x-y)^4 (x+y)^a U_u(x) U_v(y)]` (a stencil of `M`).
  - Summing over the Pieri pieces telescopes to
    `L_3(s,t,a) = sum_(j=0)^t M_2(s+t-2j,0) + M_2(s+t+2,0) - M_2(s+1,t+1)`.
    The exception is absorbed: `L_3(2,1,1) = 10 - 2 = 8`.
  - Proved: parity, support (`L_3 = 0` for `s-t > a+4`), the boundary
    `s-t = a+4` (`L_3 = 1`), and `a = 0`.  For `s+t > a+2` it reduces to
    the one-plus-label `r = 2` row sum.
  - Open: (5) `sum_j M_2(s+t-2j,0) + M_2(s+t+2,0) >= M_2(s+1,t+1)` for
    `s+t <= a+2`.
  - Main agent: the packet verifier passes, and `L_3` agrees with an
    independent Kostka/branching computation in 405 cases.  The packet's
    screen has 115,506 grouped cases, no negative value.
- **Observation (main).**  With `p = s+1`, `q = t+1`, the same telescoping
  at `r = 2` gives `L(s,t,a) = sum_(j in C(p,q)) tau_j - M(p,q;a)`, i.e.
  FM42's formula `T = sum tau_j + M` with the cross term's sign flipped
  (`D_p D_q` instead of `S_p S_q`).  So FM54 is the `D D` analogue of FM42.
- **`[q'^-, 1^-, p^+, q^+, 1^+ x a]` (venus).**
  - Decomposition `J/2 = sum_(j in C(p,q)) A_j(s,a) + X`.  Every
    `A_j >= 0` by FM53, and `B_j(lambda) = 1` iff
    `(A+C = j, min(A,C) = 0)` or `(A+C = j-2, B = 0)`.
  - `X` is an FM30 cross sum.
  - Proved: `a = 0` (`X = 0` exactly), parity, and outer support
    `|p+q-s| > a+2`.  Screen: 1,944 cases, no negative value.  `X < 0`
    already at `(2,2,2;2)`.
- **MP_2 (mars).**  No proof.
  - Screen through degree 30: 28,629 partitions, no negative value.
  - A Schur-positive `h_1^perp` Pieri induction fails: `s_(3,1,1)` gets
    `-1`.

### Corollary FM56 (two-sided row-interval bound for the `r = 1` cross kernel; from FM42 and FM54)

For all `p, q >= 1` and `a >= 0`:

```text
|M(p,q;a)| <= sum_(j in C(p,q)) tau_j(a),     M(p,q;a) = E[(x-y)^2 (x+y)^a U_p(x) U_q(y)],
```

where `C(p,q) = {|p-q|, |p-q|+2, ..., p+q}` and `tau_j = M(j,0)`.
- The `+` side (`sum tau + M >= 0`) is FM42: plus labels `p, q`.
- The `-` side (`sum tau - M >= 0`) is FM54, via
  `L(p-1,q-1,a) = sum_(j in C(p,q)) tau_j - M(p,q)`: minus labels `p, q`.
  Checked exactly in 880 cases.

**Pattern.**  Let
`M_(r-1)(u,v;a) = E[(x-y)^(2(r-1)) (x+y)^a U_u(x) U_v(y)]` and
`tau_(r-1)(j) = M_(r-1)(j,0)`.  The two-label families at level `r`,
`[p^pm, q^pm]` with fundamental labels, are exactly
`|M_(r-1)(p,q)| <= sum_(C(p,q)) tau_(r-1)(j)`.
- Neptune's (5) is the `-` side at `r = 3`.
- Jupiter's row-interval inequality (1) is the `r = 2` analogue for a
  single constituent.

- **FM-MP-C (mercury).**
  - (a) The three-term Schur formula: ACCEPTED.  Independent derivation
    of the Laplace-expansion signs `+,-,+,+,-,+` and pairwise symmetry,
    plus small cases.
  - (b) `r = 1`, `k = 2`: exact `B_p = 1` shape families (with
    determinant shifts):
    - `(p+b+d, b+d, d, d)`;
    - `(p+b+d, p+b+d, p+d, d)`;
    - `(p-2+d, c+d, c+d, d)`, with `0 <= c <= p-2`.

    The aggregate Pieri inequality over these families is unresolved.
    Shape-by-shape and K-type-by-K-type comparisons both fail, e.g. at
    `p = q = 2`, `kappa = (1,1)`.

- **FM-T1-3FA (mars): the three-factor identity and a unification.**
  - With `p = s+1`, `q = t+1`, `r = u+1`, Clebsch--Gordan on `D_p D_q`
    gives
    `L(s,t,u;a) = sum_(j in C(p,q)) F(j, r; a) - K(p,q,r;a)`, where:
    - `F(n,r;a) = (1/2) int D_1 D_r S_n (x+y)^a = <hat S_n, Sym^(r-1) W (x) W^a>`,
      which is `>= 0` by FM53;
    - `K` is the cross term `(1/2) int D_1 D_r (U_p(x)U_q(y) + U_p(y)U_q(x)) (x+y)^a`.

    Exact formulas for `F` and `K` in `M`; the packet's verifier passes.
  - **Unification (main agent).**  Venus's family
    `[r^-, 1^-, p^+, q^+, 1^+ x a]` has value `sum_j F(j,r;a) + K`, with the
    same `F` and `K`.  So the two families together are the two-sided bound
    `|K(p,q,r;a)| <= sum_(j in C(p,q)) F(j,r;a)`, where `F` is the
    FM53-positive sum.  This generalises FM56 (the case `r = 1`).  In
    general, for a polynomial GL(4)-module `E`, the families
    `[p^+, q^+ | E]` and `[p^-, q^- | E]` together say
    `|<cross_(p,q), E>| <= sum_(j in C(p,q)) <hat S_j, E>`.

- **FM-T2-PQS2 (venus): closed forms for the minus-label-3 (`s = 2`) two-sided bound.**  Here "r = 3" means minus label 3 in an `r = 1` word, not the level-3 weight `Q_3`.
  They serve both the `+` side (venus's family) and the `-` side (three
  factors with `u = 2`).
  - `A_j(2,a) = F(j,3;a) = g_a((j-3)/2) - g_a((j+5)/2) - m_a(j-1,3)` for
    `j >= 4`, via FM49 and telescoping, with
    `g_a(t) = beta_a(t)^2 - beta_a(t-1) beta_a(t+1)`.  Smaller `j` are
    given by FM53 branching.
  - `X_2 = K(p,q,3;a) = M(p,q+2) + 2M(p,q) + M(p,q-2) + M(p+2,q) + M(p-2,q) + sum_(eps,eta = +-1) M(p+eps, q+eta)`.
  - Proved: `X_2 = 0` for `p+q >= a+5`, plus parity.  Open:
    `sum_(j in C(p,q)) A_j(2,a) >= |X_2|`, i.e. both signs, for
    `a >= 1`, `p+q <= a+4`.
  - Exact checks: 182 `A_j` identities and 637 decompositions, all
    `J >= 0`.  `X_2 < 0` already at `(2,3,1)` and `(2,2,2)`.

- **FM-T2-R2PM3 (jupiter): the `r = 2` one-plus-label family as a window inequality.**
  - With `a = d+2n`, `x = n-q` and FM54's `A, B, C`, the family is exactly
    `D_(n+1) - D_x >= |A(b_(x+1) - b_(x-4)) + B(b_x - b_(x-3)) + C(b_(x-1) - b_(x-2))|`.
    Here `c_j` is the third difference of the binomial row, and
    `D_j = c_j^2 - c_(j-1) c_(j+1)`.
  - Parity and off-support cases: proved.
  - Screen: 597,861 active triples (`a <= 240`), minimum slack 2.
  - Kills:
    - `R >= Phi` fails (at `(a,d,q) = (5,5,1)`, `Phi = 7 > R = 6`);
    - the gap `R + Lambda` is not unimodal in `x`.
  - Main agent: every increment `D_(y+1) - D_y` below the centre is
    nonnegative (`a <= 120`).  But `D_(x+1) - D_x >= |Lambda|` fails in
    76,421 of 77,531 cases, so the whole telescoped window is needed.
  - Symbolic: `(D_(x+1) - D_x)/b_x^2` equals a positive factor times an
    explicit quartic `P_a(x)`.

- **FM-T1-3F3 (saturn): recursion in the suffix length.**
  - For `M_a = Sym^s (x) Sym^t (x) Sym^u (x) W^a`, applying the `W`-Pieri
    rule twice gives
    `L(a+2) - L(a) = -2m_00 + 4m_11 + m_20 - m_31 - 2m_22 + m_40` (at `M_a`).
  - With the base cases `a = 0, 1` (FM55), the whole three-factor family
    reduces to the increment inequality (4), `Delta >= 0`.  (4) is OPEN;
    this is a reduction, not a proof.
  - Screen: 4,781 increments (labels `<= 14`, `a <= 20`), none negative.
  - Kill: decomposition into FM54 prefix intervals fails, since
    multiplicities are not monotone along diagonals (`Sym^12` cubed:
    `m_10,0 = 66 < m_11,1 = 143`).
  - Note (main): (4) is an FM3-type statement for the virtual weight
    `Q_2 (x) (W (x) W - 1) = Q_2 (x) (V_(2,0) + V_(1,1))`.  In the multigraph
    model each count `I` increases when two leaves are appended (a
    leaf--leaf edge at the right end creates no 3-nesting).  But `L` is a
    signed combination, so this alone does not give (4).
- **Conjecture MON_2 (main agent, exact screen).**  For every nonempty
  composition `kappa`, `L(kappa (+) (1,1)) >= L(kappa)` at `r = 2`: adding two
  unit factors never decreases the H-only value.
  - Screen: all 272 partitions of size `<= 12`, with exactly one failure,
    `kappa = ()` (`L = 5` versus `L(1,1) = 3`).
  - MON_2 reduces FM3 at `r = 2`, H-only, to products with at most one
    unit factor.  Saturn's (4) is its three-factor case.
- **FM-T1-3FB (mars): minus side at `r = 3`.**
  - The same closed forms as venus, derived independently: `F(j,3;a)` by
    FM49 plus telescoping, or as a six-term `M` stencil (the edge
    indicators matter at `j = 0, 1`); `K(p,q,3;a)` as the ten-term stencil.
  - The telescoped `sum_(j in C(p,q)) F` leaves at most eight `g_a`
    boundary terms.
  - Verifier passes (143 `F` checks, 468 `K` and fusion checks).
  - The sign `K <= sum F` is open for `a >= 2`, `s, t >= 2`.

### Level-mixing identity for the `u = 2` slice (main agent)

In the doubled ring, `xy = ((x+y)^2 - (x-y)^2)/4`, so
`H_3 = (x+y)^2 - xy - 2 = (3u^2 + v^2)/4 - 2`, with `u = x+y` and
`v = x-y`.  Hence, exactly,

```text
L(s,t,2;a) = (3/4) L_2(s,t;a+2) - 2 L_2(s,t;a) + (1/4) L_3(s,t;a),
```

where `L_2` is the FM54 two-power value at `r = 2` and `L_3` is the
two-power value at `r = 3` (neptune's family).  Checked exactly in 405
cases.

**Reduction (CONDITIONAL: both (i) and (ii) are open).**  For `a >= 2`, the `u = 2` slice follows from:
- (i) `L_3(s,t;a) >= 0`, which is neptune's (5);
- (ii) the growth bound `3 L_2(s,t;a+2) >= 8 L_2(s,t;a)` for `s >= t >= 2`,
  `a >= 2`.  It FAILS at `(s,t,a) = (0,0,2)`: ratio 2 (FM-T1-GR, saturn).
  In an exact screen of 1,450 cases (`s <= 16`, `2 <= a <= 22`) that is the
  only failure, and none occur for `t >= 1`.

Evidence for (ii): the minimum of `L_2(a+2)/L_2(a)` over `s >= t >= 2`
is 3.0 at `a = 1` and 3.75 at `a = 2`, and it increases roughly like
`a/2 + 2` (labels `<= 14`, `a <= 14`).  It fails only at `a = 0`
(ratio `4/3`), which FM55 covers.  With (i), the slack is still at least
`4/3` at `a = 0` over the screen.

More generally, every extra `H` factor is a polynomial in `(u, v^2)`, so
every H-only product is a signed combination of two-factor values at
shifted `(r, a)`.  This is the note's earlier `(U,V)` and `E_(m,r)`
machinery.
- **FM-T1-3F4 (saturn): the increment as a `K`-coefficient inequality.**
  - Identity: `Delta(M) = -2c_00 - 5c_02 + 3c_04 + c_06 + 6c_11 + 4c_13 - 2c_15 - 6c_22 - c_24 + 2c_33`.
    The target `Delta(M) >= 0` is OPEN and intended only for nonempty
    products; it is false for the empty product, see below.
    where `c_pq` are the SU(2) x SU(2) coefficients of `Res_K M`.  It
    comes from `(x-y)^4 ((x+y)^2 - 1)`.
  - Kills:
    - the empty product (`Delta(1) = -2`);
    - per-Schur positivity: `Delta(S_(4,2,2)) = -2`, while
      `Delta(h_4 h_2 h_2) = 3`.
  - Screen: 209 products of 1--4 factors, minimum 0.

### FM-CHK16 (luna_max_mercury): FM55, FM56, the `r = 3` reduction

- **FM55: mathematics ACCEPTED; two note repairs applied.**
  - `B, D, C` are now defined, and the open range is stated as `a >= 2`.
  - Mercury re-derived the invariant-ring input from Weyl's symplectic
    first and second main theorems: the ring is free for at most five
    slots, and for six slots it has the single Pfaffian relation.
  - The Pfaffian is a nonzerodivisor, since the pairing ring is a domain.
  - The complete case checks for `a = 0, 1` and the parity, support and
    unit-factor reductions all check.
- **FM56: ACCEPTED.**  Clebsch--Gordan in the doubled ring gives
  `S_p S_q` and `D_p D_q` as `sum S_j +- cross`, with small checks.
- **`r = 3` reduction (neptune): ACCEPTED.**  Checked `Phi_3`, the
  telescoped `L_3`, the support and boundary values, and small cases.
  Interior (5) remains open.
- **FM-T2-PQS3 (venus): further telescoping, no sign.**
  - `m_a(j-1,3) = Phi_a(t) - Phi_a(t+1)`, where
    `Phi_a(t) = beta_a(t+3) beta_a(t) - beta_a(t+4) beta_a(t-1)`.
  - So `sum_(j>=4) A_j(2,a)` is a four-endpoint `g_a` window plus two `Phi_a`
    endpoints.  All terms are shifted rows of one binomial sequence
    `B_r = binom(a+2, r)`.
  - Screen: 380 active cases, minimum 2.  The sign comparison is open.
- **FM-T1-TR3r2 (neptune): partial.**
  - Identity: the margin of (5) equals the original integral
    `(1/2) E[(x-y)^6 (x+y)^a H_(s+1) H_(t+1)]`.
  - Proved:
    - the square case `s = t` (integrand `>= 0`);
    - the strips `t = 0, 1` with `s <= 3`, from FM18's `S^(2m) T^(2r)`
      moment formula and `E[S^(2m) T^(2r) H_3]`.
  - Open: `s > t` in general.  The first open cases are `(4,0,2)`
    (margin 7) and `(3,2,3)` (margin 15).
- **FM-T2-R2PM4 (jupiter): no progress; the line is paused.**
  - The two signs of the row-interval inequality (1) are exactly the word
    families `[(P+1)^-, 1^-, 1^-, 1^-, Q^+, 1^+ x a]` (`R + Lambda`) and
    `[Q^-, 1^-, 1^-, 1^-, (P+1)^+, 1^+ x a]` (`R - Lambda`), with the
    factor 2 from FM6.
  - `Lambda(P,Q) = Phi(P-Q,Q) - Phi(P-Q,Q-1)`.
  - The `Q = 1` plus sign equals `Lambda(P+1,0) + Phi(P-1,1) >= 0`, by FM19
    and FM54 (already inside FM19's family).
  - Screen: 1,496 active triples (`a <= 30`), all nonnegative.

## Notation guide for the FM47--FM56 section (added after FM-AUDIT1)

Several symbols are reused locally.  Each entry's own definitions apply.
- `Phi`: the rectangle ratio in FM47 (`Phi_rect`); the two-row partial
  sums `Phi(d,t,a)` in TR/FM54; `Phi_a(t)` in FM-T2-PQS3.
- `M`: the FM18 normaliser `M(m,r)` in FM47; the FM30 kernel `M(u,v;a)`
  from FM54 onwards; the moment matrices `M^k_ij` in FM-L123m; the module
  `M` in FM55 and FM-T1-3F.
- `beta`: the ratio `b/N` in FM47; `beta_a(t)` in FM49; `beta` in FM54;
  `beta_x = binom(N, N/2+x)` in FM30; `|c-d|` in FM50.
- `Delta`: `Delta_(p,q)(a)` in FM49; the SO(5) functional `Delta(M)` in
  FM50; the increment `Delta(M)` in FM-T1-3F3/4.
- `Q`: the SO(5) module `Q` in FM50; the Sp(4) weights `Q_r`; the
  polynomials `Q_1, Q_2` in FM49.
- `A, B, C`: the window coefficients in FM54; part differences in the
  three-term formula; indicator sums in FM55.

## Appendix to FM54 (verbatim from FM-T1-TR, luna_max_saturn; inserted after FM-AUDIT1)

Repairs from FM-CHK-TR apply.
- `R_0` has discriminant 144 at `d = 0, 1`, where `R_0 = 12(n+1)(n+2)`.
- The `n = 0` endpoint: `delta > 0` forces `d >= 4`, and
  `F(1) - F(0) = (d^3 - 2d^2 - 5d + 30)/(2(d+3)) > 0`.
- Section 5 below treats `n = 0` separately in any case.

#### 2. The tail-window formula

Write bⱼ = C(a,j), with bⱼ = 0 outside 0 ≤ j ≤ a. If a and d have different parity, every Weyl-binomial term vanishes, so Φ(d,t,a) = 0. Suppose they have the same parity and d ≤ a; put n = (a−d)/2, so a = d+2n.

For G(d,s,a) = Σ₍q≥s₎ Λ(q+d,q;a), pairing each Weyl term indexed by (u,v) with the term indexed by (−v,−u) telescopes its q-sum to a window of binomials. Substituting the three targets (2,1), (3,2), (4,1), with weights 5, −3, 1, gives

G(d,s,a) = A Σᵢ₌₀⁴ b₍d+n−1+s+i₎
      + B Σᵢ₌₀² b₍d+n+s+i₎ + C b₍d+n+1+s₎,

where

A = −3(bₙ−bₙ₋₁) + bₙ₊₁−bₙ₋₂,

B = 5(bₙ−bₙ₋₁) + bₙ₋₃−bₙ₊₂,

C = 5(bₙ₋₂−bₙ₊₁) − 3(bₙ₋₃−bₙ₊₂).

This is the finite window identity recorded in the note’s telescoping paragraph at [line 3989](/home/yang/q3adjoint/ginibre_q3/SU2_FUNDAMENTAL_MINUS_REDUCTION_2026_09_26.md:3989).

By binomial symmetry b₍a−j₎ = bⱼ, set x = n+1−s. Normalize A, B, C by bₙ and call the resulting coefficients A₀, B₀, C₀. Define

α = A₀, β = B₀−3A₀, γ = A₀−B₀+C₀,

and

F(x) = α C(a+4,x) + β C(a+2,x−1) + γ C(a,x−2).

Expanding
A₀(1+z+z²+z³+z⁴) + B₀(z+z²+z³) + C₀z²
in the basis (1+z)⁴, z(1+z)², z² gives the exact identity

G(d,s,a) = bₙ F(n+1−s).                                      (1)

The ratios bₙ₋ⱼ/bₙ and bₙ₊ⱼ/bₙ, for j = 1,2,3 as needed, are respectively

n/(d+n+1), n(n−1)/((d+n+1)(d+n+2)), 
n(n−1)(n−2)/((d+n+1)(d+n+2)(d+n+3)),

and

(d+n)/(n+1), (d+n)(d+n−1)/((n+1)(n+2)).

Substitution gives

α = (d+1)δ / ((n+1)(d+n+1)(d+n+2)), δ = d²−d−6n−6.       (2)

A second simplification of the window formula gives

G(d,0,a) = bₙ² (d+1)(d+2n+1)(d+2n+2) R₀
 / ((n+1)²(n+2)(d+n+1)²(d+n+2)²(d+n+3)),                 (3)

where

R₀ = d⁴+2d³−4d²n−5d²+4dn+2d+12n²+36n+24.

In particular G(d,0,a) > 0. Indeed, as a quadratic in n, R₀ has leading coefficient 12 and discriminant
−16(2d⁴+8d³+2d²−12d−9). For d ≥ 2, the polynomial in parentheses is positive (its value at 2 is 71 and it increases thereafter), so the discriminant is negative. For d = 0,1, both remaining coefficients are positive.

#### 3. The endpoint comparison

For n ≥ 1, we have G(d,n+1,a) = bₙF(0) = bₙα and G(d,0,a) = bₙF(n+1). I claim

F(n+1) ≥ F(0).                                                (4)

If δ ≤ 0, this follows at once from (2) and (3). If δ > 0, then

F(n+1)/α =
bₙ (d+2n+1)(d+2n+2)R₀
/ ((n+1)(n+2)(d+n+1)(d+n+2)(d+n+3)δ).

Use the identity

R₀−δ² = 2(d+2n+1)(δ+d²+2d).

It implies R₀/δ ≥ 2(d+2n+1). For n ≥ 2, the three factors in the denominator satisfy

(d+n+1)(d+n+2)(d+n+3) ≤ (d+2n+1)²(d+2n+2).

Also bₙ = C(d+2n,n) ≥ C(2n,n) ≥ (n+1)(n+2)/2. Thus F(n+1)/α ≥ 1.

For n = 1, δ > 0 forces d ≥ 5; direct substitution reduces the ratio to R₀/(6δ). Here
R₀−6δ = d²(d²−15)+2d³+12d+144 > 0,
so (4) holds in this case too.

#### 4. The maximum of F on the required indices

For 1 ≤ x ≤ n, define

P(x) = α(a+4)(a+3)(a+2)(a+1)
 + β(x+1)(a+4−x)(a+2)(a+1)
 + γx(x+1)(a+4−x)(a+3−x).

The binomial-ratio formulas above give the exact identities

F(x+1)−F(x) =
(a+3−2x) C(a+4,x) P(x)
 / ((x+1)(a+4)(a+3)(a+2)(a+1))

and

P(x) =
(d+1)(d+2n+1)(d+2n+2)(x−n−1)(x−d−n−2) Ψ(x)
 / ((n+1)(n+2)(d+n+1)(d+n+2)(d+n+3)).

To specify the remaining quadratic, put y = n+1−x and set

Ψ(x) = 4εy² + (4d³−24dn−52d−24n−48)y
 + d⁴−2d³−12d²n−17d²+36dn+66d+60n²+228n+216,

where ε = d²−d−6n−12. These identities follow by expanding the displayed P(x) after substituting the five binomial ratios.

For x ≤ n, all factors outside Ψ(x) in the difference formula are positive: in particular a+3−2x ≥ d+3, and the two factors x−n−1 and x−d−n−2 are both negative. Therefore the sign of F(x+1)−F(x) is the sign of Ψ(x).

At x = n,

Ψ(n) = d⁴+2d³−(12n+13)d²+(12n+10)d+60n²+180n+120 > 0.    (5)

For d = 0,1,2,3,4, this is a quadratic in n with positive coefficients (its constant terms are 120, 120, 120, 168, 336). For d ≥ 5, its discriminant as a quadratic in n is
−48(2d⁴+16d³+22d²−40d−75) < 0,
and its leading coefficient is 60. This proves (5).

As a function of y, Ψ has leading coefficient 4ε. If ε > 0, its linear coefficient is positive as well:
4d³−24dn−52d−24n−48
= 4d(d²−6n−13)−24(n+2) > 0.
Indeed ε > 0 gives d²−6n−13 > d−1 and n+2 < (d²−d)/6. Hence Ψ increases for y ≥ 1 and is positive there by (5). If ε ≤ 0, Ψ is concave or linear in y. Since Ψ(n)>0, its sign as y runs from 1 upward can change only from positive to nonpositive. Reversing the order, the signs of F(x+1)−F(x), for x = 1,…,n, can change only from nonpositive to nonnegative. Thus on x = 1,…,n+1, F has no interior maximum: its maximum is at x = 1 or x = n+1.

It remains to account for x = 0. From (2) and the coefficient formulas,

F(1)−F(0) =
(d+1)H / ((n+2)(d+n+1)(d+n+3)),

where H = d³+2d²n−2d²−8dn−5d−12n²−6n+30. Also

Ψ(1) = (d+2n+3)(H+3(8−δ)).

If F(1)>F(0), then H>0. When δ ≤ 8, this gives Ψ(1)>0. When δ > 8, we have d ≥ 5 and

H+3(8−δ) = (δ(d²+2d−18−δ)+144)/3 > 0,

because δ ≤ d²−d−6 implies d²+2d−18−δ ≥ 3d−12 > 0. Thus Ψ(1)>0 in either case. If ε ≤ 0, concavity and the positivity of Ψ at x = 1 and x = n imply Ψ(x)>0 for 1 ≤ x ≤ n. If ε > 0, the preceding increasing-in-y argument gives the same conclusion. Hence F(x) increases from x = 1 to x = n+1, and F(0)<F(n+1).

Combining the cases with (4), we have proved

F(x) ≤ F(n+1) for every integer 0 ≤ x ≤ n.                    (6)

#### 5. Conjecture TR and boundary cases

For n ≥ 1, equation (1) and (6) give G(d,s,a) ≤ G(d,0,a) for 1 ≤ s ≤ n+1. For s ≥ n+2, equation (1) gives G(d,s,a)=0, while (3) gives G(d,0,a)>0. Therefore

Φ(d,t,a) = G(d,0,a)−G(d,t+1,a) ≥ 0.

When n = 0, equation (1) gives G(d,s,a)=0 for s ≥ 2, while (3) simplifies to

G(d,0,a) = (d⁴+2d³−5d²+2d+24)/(2(d+2)(d+3)) > 0.

For t = 0, Φ(d,0,a)=Λ(d,0;a)≥0 is exactly the FM19 one-row case with s=d=a. Its hypotheses match; FM19 gives the small values 5, 3, 2, 2 for a=0,1,2,3 and its stated positive formula for a≥4. For t ≥ 1, Φ(d,t,a)=G(d,0,a)>0.

If d > a, parity mismatch gives zero unless d−a is even. In the matching-parity cases d ≥ a+2. The Weyl coefficient is zero when its displacement has ℓ₁-length greater than a. For q ≥ 1, the minimum such length to the three target orbits is at least d+2q−2. For q=0 the three lengths are d, d, and |d−2|. Consequently, if d=a+2, only the (2,0) target at q=0 contributes; FM19’s boundary value is 1, so Φ(d,t,a)=1 for every t. If d ≥ a+4, all terms vanish, so Φ(d,t,a)=0. This completes the remaining support cases.

The smallest values, computed from the Weyl sum, are:

| (a,d,n) | Φ at t=0 | Φ at t=1 | later values |
|---|---:|---:|---:|
| (0,0,0) | 5 | 2 | 2 |
| (1,1,0) | 3 | 1 | 1 |
| (2,2,0) | 2 | 1 | 1 |
| (3,3,0) | 2 | 2 | 2 |
| (2,0,1) | 3 | 4 | 2 from t=2 onward |

For (a,d)=(2,0), the successive Λ terms are 3, 1, −2, followed by zeros. The outer boundaries give Φ=1 when d=a+2, Φ=0 when d≥a+4, and Φ=0 for parity mismatch.


## Appendix to FM50 (verbatim from FM-L123k, luna_max_neptune; inserted after FM-AUDIT1)

**2. Weyl-kernel lemma.** Define
\[
\kappa(t)=
\begin{cases}
2,&t=0,\\
-1,&|t|=1,\\
0,&|t|\ge2.
\end{cases}
\]
For irreducibles \(\lambda=[x,y]\), \(\nu=[r,s]\), use coordinates
\[
U=x+y,\quad V=x-y,\qquad
A=r+s+2,\quad B=r-s+1.
\]
Then the coefficient \(K(\lambda,\nu)=\langle\lambda\nu,Q\rangle\) is obtained from the Weyl alternant by applying the finite weight stencil
\[
\omega(0,0)=4,\quad
\omega(\pm2,0)=\omega(0,\pm2)=-2,\quad
\omega(\pm2,\pm2)=1,
\]
with all other \(\omega\) zero.

Indeed, for \(X=z_1+z_1^{-1}\), \(Y=z_2+z_2^{-1}\), the SO(5) characters are
\[
\chi_{[1,0]}=1+X+Y,\quad
\chi_{[1,1]}=2+X+Y+XY,\quad
\chi_{[2,0]}=X^2+Y^2+XY+X+Y-2.
\]
Thus \(\chi_Q=(X-Y)^2\), which has exactly the displayed weights in \((U,V)\)-coordinates. The Weyl character formula says that \(K(\lambda,\nu)\) is the coefficient of the monomial of exponent \(\lambda+\rho\) in \(\chi_Q A_{\nu+\rho}\), where \(A_{\nu+\rho}\) is the alternating Weyl sum. This is Fulton–Harris, *Representation Theory: A First Course*, §24.1; its finite-dimensional semisimple Lie algebra hypotheses apply to \(\mathfrak{so}_5(\mathbb C)\), of type B₂. [Section reference](https://www-fourier.ujf-grenoble.fr/~panchish/ETE%20LAMA%202018-AP/lecturesZETAS2018/Fulton%20Harris_Representation%20theory%20first%20course.pdf)

**3. Exact formula.** Put
\[
\alpha=|a-b|,\quad m=\min(a,b),\qquad
\beta=|c-d|,\quad n=\min(c,d),
\]
and for integer intervals define
\[
\operatorname{ov}([r,s],[u,v])
=\max(0,\min(s,v)-\max(r,u)+1).
\]
Then:

- If \(\alpha\equiv\beta\pmod 2\), put \(h=(\alpha-\beta)/2\). Then
\[
\begin{aligned}
\Delta(H_aH_bH_cH_d)={}&
\operatorname{ov}([h,h+m+1],[0,n+1])\\
&+\operatorname{ov}([h,h+m],[0,n])\\
&+\mathbf1_{\{h=0,\ m=n\}}
+\mathbf1_{\{\alpha=\beta=0\}}\bigl(1+\mathbf1_{\{m=n\}}\bigr).
\end{aligned}
\]

- If \(\alpha\not\equiv\beta\pmod 2\), put \(h=(\alpha-\beta-1)/2\). Then
\[
\begin{aligned}
\Delta(H_aH_bH_cH_d)={}&
\operatorname{ov}([h+1,h+m+1],[0,n+1])\\
&+\operatorname{ov}([h,h+m+1],[0,n]).
\end{aligned}
\]

**Proof.** The pair decomposition recorded in FM44 is
\[
H_aH_b=\bigoplus_{J=0}^{m}\ \bigoplus_{K=0}^{m-J}
[a+b-2J-K,K].
\]
Writing \(i=m-J,\ j=i-K\), its support in \((U,V)\)-coordinates is
\[
(U,V)=(\alpha+2i,\alpha+2j),\qquad 0\le j\le i\le m.
\]
The other pair has the same form with \(\beta,n,k,l\). Since \(\rho=(3/2,1/2)\), the shifted points are
\[
(\alpha+2i+2,\alpha+2j+1),\qquad
(\beta+2k+2,\beta+2l+1).
\]

The Weyl group acts by signed permutations in these coordinates. The weight stencil in Statement 2 has coefficients \(\kappa(u)\kappa(v)\) at \((2u,2v)\). Matching the shifted points against the signed Weyl orbit gives these pair contributions:

- For equal parity, with \(h=(\alpha-\beta)/2\), the identity Weyl element contributes
\[
\kappa(h+i-k)\kappa(h+j-l).
\]
The only other possible element is the reflection of the second coordinate. It contributes
\[
\mathbf1_{\{\alpha=\beta=0,\ j=l=0\}}\kappa(i-k).
\]
All other signed or swapped images either have the wrong parity or differ by at least \(3\) in a coordinate, outside the stencil’s support.

- For opposite parity, the only possible Weyl element is the positive coordinate swap. With \(h=(\alpha-\beta-1)/2\), its contribution is
\[
-\kappa(h+1+i-l)\kappa(h+j-k).
\]
The other swapped images have a coordinate difference at least \(3\); the non-swapped images have the wrong parity.

To sum these terms, let \(e_t\) be the standard basis of finitely supported sequences on \(\mathbb Z\), and put \(b_t=e_t-e_{t+1}\). Then \(\langle b_r,b_s\rangle=\kappa(r-s)\). Define
\[
B_m=\sum_{0\le j\le i\le m}b_i\otimes b_j
=\sum_{r=0}^{m+1}e_r\otimes e_r
-\sum_{r=0}^{m}e_r\otimes e_{r+1}
-e_{m+1}\otimes e_0.
\]

In the equal-parity case, the identity contributions sum to \(\langle\tau_hB_m,B_n\rangle\). The diagonal and adjacent-diagonal overlaps give the first two interval terms in the formula; the corner terms coincide exactly when \(h=0,m=n\). The reflection contributions sum to
\[
\sum_{i=0}^m\sum_{k=0}^n\kappa(i-k)
=\left\langle e_0-e_{m+1},e_0-e_{n+1}\right\rangle
=1+\mathbf1_{\{m=n\}}.
\]
This gives the equal-parity formula.

In the opposite-parity case, the signed sum is
\[
-\left\langle(\tau_{h+1}\otimes\tau_h)B_m,\ B_n^{\mathsf T}\right\rangle.
\]
Expanding the displayed expression for \(B_m\), the only overlaps are the main-diagonal and adjacent-diagonal segments. Their signed counts are the negatives of the two overlap terms in the formula, proving the opposite-parity expression.

Every term in both formulas is nonnegative, proving the theorem. ∎


### FM-AUDIT1 (luna_max_venus): consistency audit of the FM47--FM56 section

The repairs listed in the audit have been applied:
- the FM47 raw-loss sentence;
- stale checker and open-range text for FM47, FM49, FM54, TR and 2L11;
- the three-term section's final `>= 0` is now marked OPEN;
- the increment `Delta >= 0` is now marked OPEN (false for the empty
  product);
- the level-mixing reduction is marked CONDITIONAL;
- definitions added: FM48 `J1, J2, R`; FM49 parity and `Q_1, Q_2`, and
  the verifier filename; FM50 `beta, n` and the group convention; the
  FM51 lattice parity; FM55 `R`, `B(k-(1,1,1))` and the check bounds;
  FM53's `j = 0` edge; the local notation for FM47;
- the notation guide above;
- verbatim appendices for FM54 and FM50, so that both proofs can be
  reproduced;
- the "r = 3" wording in PQS2/3.

**Final status (audit).**
- Proved and independently checked: FM47, FM49, FM51, FM53, FM54,
  FM55 (`a <= 1` and the other listed slices), FM56, and the three-term
  Schur formula.
- FM50: proved, and independently checked for all labels (FM-CHK17).
- FM48: proved in the support range (main agent only).
- FM52: proved (a consequence of FM49).
- Conjectures: MP, MP_2, MP_r, MON_2.
- Reductions and screens only: FM-L123l/m, FM-MP-B, FM-T1-3F3/4,
  FM-T1-3FB, FM-T2-PQS2/3, and the `r = 3` interior.

- **Kill K26 (FM-T1-3FC, mars): no constant-coefficient certificate for
  the `u = 2` slice.**
  - Tested 106 shifted proved blocks: 39 FM54 values, 15 FM53 values,
    35 FM42 values, 12 `Phi` values and 5 pure rows.
  - An exact Farkas certificate (16 target points in the open region)
    shows that no nonnegative constant-coefficient combination of these
    blocks equals `T = L(s,t,2;a)` with zero remainder.
  - This does not exclude parameter-dependent coefficients or a positive
    remainder.  A fitted identity `T = X + 407 F(s+t+2,0;a+2)` fails at
    `(2,2,2)`.

### FM-CHK17 (luna_max_mercury): FM50, all labels -- ACCEPTED

Mercury checked every step of the appendix proof:
- the kernel stencil;
- the pair supports from FM44's all-label SO(5) product formula;
- the classification of the Weyl images: excluded images are at distance
  `>= 4`, and only the identity and the `V`-reflection (the latter when
  `alpha = beta = 0`) or the swap `(V,U)` contribute;
- the `B_m` identity and the overlap formulas, including the boundaries
  `m = 0`, `n = 0`, `h = -1, 0`;
- nonnegativity, with small cases.

There is also a supplemental exact enumeration (`alpha, beta <= 7`,
`m, n <= 5`).  FM50 now has an independent all-label check.
- **FM-T1-GR2 (saturn): growth bound reduced to a termwise estimate.**
  - For `s >= t >= 2`, (G) follows from
    (TG) `Phi(d,j,a+2) >= 3 Phi(d,j,a)` for all `d, j >= 0` with
    `d+2j >= 4`, together with FM54.
  - Proved for the support cases (`d >= a+2`, via FM19 and FM54's `n = 0`
    formula).
  - The interior is the explicit window inequality
    `binom(a+2,n+1)[F_(d,n+1)(n+2) - F_(d,n+1)(n+1-j)] >= 3 binom(a,n)[F_(d,n)(n+1) - F_(d,n)(n-j)]`.
  - Screen: (G) in 8,004 cases (minimum ratio `15/4`), and (TG) in
    82,708 summands; no failures.
- **FM-T1-TR3r3 (neptune): the FM54 machinery at `r = 3`.**
  - Window formula `G_3 = K_7 W_7 + K_5 W_5 + K_3 W_3 + K_1 W_1`.
  - Four-row binomial form
    `F_(d,n)(x) = alpha_3 binom(a+6,x) + beta_3 binom(a+4,x-1) + gamma_3 binom(a+2,x-2) + delta_3 binom(a,x-3)`,
    with `Phi_3 = b_n (F(n+2) - F(n+1-t))`.
  - Exact increment factorisation through a quartic `Psi_3(y)`, with
    explicit `P_0..P_3`.
  - The `n = -1` boundary is proved:
    `Phi_3(a+2,0,a) = (a^2 - 9a + 28)/2 > 0` and
    `Phi_3(a+2,t,a) = (a^2 - 7a + 18)/2 > 0` for `t >= 1`.
  - Main agent: unlike FM54, the increments change sign TWICE in `x`,
    in the pattern up, down, up (e.g. `d = 0, n = 8`:
    `+ + + + - - - + + +`).  So the proof must show that the endpoint
    `F(n+2)` beats an interior local maximum, not only `F(0)`.
  - Exact screen: `max_x F = F(n+2)` for all `d, n <= 40` (1,681 pairs),
    except `(d,n) = (1,0)`, the grouped exception.  In 354 pairs the
    sign pattern has two changes.
- **FM-T1-GR3 (saturn): (TG) for `n = 0, 1`; `n >= 2` reduced.**
  - Exact margins `H_(d,n,j) = Phi(d,j,a+2) - 3 Phi(d,j,a)` for `n = 0`
    (three rows) and `n = 1` (four rows).  Each is a polynomial with all
    coefficients positive after the shift `d = d_0 + x`, so (TG) holds
    for `n <= 1`.
  - For `n >= 2`, (TG) is equivalent to the normalised inequality
    `kappa^2 C_(d,n+1) - 3 C_(d,n) >= rho_(d,n,j) (kappa^2 theta R_(d,n+1,j) - 3 R_(d,n,j))`,
    plus a tail case `j = n+1` and the axis comparison for
    `j >= n+2`.
  - Screen: 72,520 cases (`d, n <= 50`), no failure.  (TG) is TIGHT: the
    ratio equals 3 at `(d,n,j,a) = (0,1,3,2)`.  Note (G) only needs
    `8/3` after summing.
  - Line paused.  Open: (TG) for `n >= 2`, and therefore (G).
- **FM-T1-TR3r4 (neptune): sign-shape theorem and two proved slices at `r = 3`.**
  - *Identity.*  `P_2 - P_1 = 2(d+1)^2 P_3`, hence
    `Psi_3(y) = P_0 + 8 P_1 z + 16 P_3 z^2` with `z = y(y+d+1)`.  This is a
    quadratic in the increasing variable `z`, so the increments of `F`
    change sign at most twice.  The maximum of `F` on `[0, n+2]` is
    therefore at an endpoint or at the single interior local maximum.
    Main agent: verified both identities symbolically.
  - *Exact gap formula.*
    `F(n+2) - F(n+2-y) = binom(m,n+2) G_(d,n)(y)`, with `m = d+2n+6` and
    `G` explicit in a rational function `R(t)`.  The endpoint value is
    `F(n+2) = binom(d+2n,n) (d+1)(d+2n+1)(d+2n+2) R_8(d,n) / D(d,n)`.
  - *Proved.*  The slices `n = 0` and `n = 1` for every `d` and `t`
    (explicit polynomial gaps; the only negative value is the grouped
    exception `Phi_3(1,1,1) = -2`).  Hence `L_3(s,t,a) >= 0` on the
    slices `h + t <= 1` and on `h = -1, t <= 2`, where
    `h = (a-s-t)/2`.
  - *Open.*  `R_8(d,n) >= 0` and `G_(d,n)(y) >= 0` at `y = n+2` and at the
    interior-maximum `y`, for `n >= 2`.
  - Screen: `d, n <= 100` (10,201 pairs); the interior maximum is at most
    `11/18` of the endpoint.
- **FM-T1-TR3r5 (neptune): partial `R_8` certificate.**
  - `R_8(d,n) > 0` on the cone `d >= 3n`, `n >= 2`: after `N = n-2`,
    `U = d - 3n`, every coefficient is positive.  Also for all `d` at
    `n = 2, 3`, by the cone plus direct rows.
  - Open:
    - `R_8 >= 0` for `n >= 4`, `d < 3n`;
    - the gap comparisons `G_(d,n)(y) >= 0` at the endpoint and at the
      interior maximum, for `n >= 2`.
  - The `r = 3` two-power family therefore remains open on the grouped
    region `h >= 0, h+t >= 2`, or `h = -1, t >= 3`.

### `r = 3` two-power family: proof plan after the sign-shape theorem (main agent, exact data)

With `F = F_(d,n)` as in FM-T1-TR3r3/4 and `a = d+2n`, the family
(`L_3 >= 0`) follows from four pieces:

- (a) `Psi_3(1) > 0`, i.e. the last increment is positive.  This is a
  polynomial inequality in `(d,n)`.  Data: `F(n+1) < F(n+2)` in every
  case with `d, n <= 80`.
- (b) At most two sign changes of the increments.  **Proved**
  (FM-T1-TR3r4).  With (a), the maximum of `F` on `[0, n+2]` is at
  `n+2`, at `0`, or at the unique interior local maximum `x_1`.
- (c) `F(n+2) >= F(0)`, except at `(d,n) = (1,0)` (the grouped exception).
  This is a closed-form endpoint comparison of FM54 type.  Data: no
  other failure for `d, n <= 80`; equality at `(0,0)`.
- (d) `F(n+2) >= F(x_1)`.  Data: the worst ratio is `11/18`, at
  `(d,n,x) = (1,3,1)`.
  - The interior maximum sits at `y_1 = n+2-x_1 ~ 2 sqrt(n)`: 6 at
    `n = 8`, 13 at `n = 32`, about 16 at `n = 64`.
  - There the binomial factor is about `5e-3` and the `R` factor about
    30--90, so `F(x_1)/F(n+2)` tends to about 0.28 for fixed `d` as
    `n -> infinity`.  It is small for large `d`.
  - So an FM42-style certificate (a scaling-limit bound for large `n`
    plus exact per-`n` polynomial checks in `d` for small `n`) is
    feasible.

**Refinement of piece (d) (main agent, exact and floating screens).**
- *Binomial bound* (neptune, proved):
  `C(m,n+2-y)/C(m,n+2) = prod_(i=1)^y (n+3-i)/(d+n+4+i) <= exp(-kappa w)`,
  with `w = y(d+y+2)` and `kappa = 2/(d+2n+7)`.
- `F(n+2-y)/F(n+2) = [binomial ratio] * P(w)`, where
  `P(w) = R(t_0-w)/R(t_0)` is an explicit cubic in `w`.
- Below the first sign change `y_-` of `Psi_3` (smaller root `z_-`), the
  increments are positive, so `F(n+2-y) <= F(n+2)` by monotonicity.  In
  885 of the screened pairs with `n >= 4` there is no sign change at
  all.
- For `y >= y_-` and `n >= 4`, `exp(-kappa w) P(w) <= 1` holds in every
  screened case.  Maximum 0.931 at `(d,n,y) = (1,4,5)`; per `n`: 0.867
  (`n = 5`), 0.75 (`n = 6`), about 0.5 (`n = 16`), about 0.335
  (`n = 400`).
- So piece (d), and piece (c) for `n >= 4`, reduce to a one-variable
  estimate with at least 7% slack:
  `exp(-kappa w) P(w) <= 1` for `w >= w_-(d,n)`.
- For `n = 2, 3`, the Gaussian bound is too weak in five pairs.  Use the
  exact binomial ratios instead: finitely many rational inequalities in
  `d` for each `(n, y)`.
- A Taylor-4 plus AM--GM certificate for `P(w) <= e^(kappa w)` on all
  `w >= 0` fails at `d = 0`, since `kappa w` reaches about `n`.  A range
  split in `w` is needed (FM42 style).

**`r = 3` pieces (a), (c), and `R_8` (FM-T1-TR3-AC, luna_max_saturn; main-agent additions).**
- **(a) PROVED: `Psi_3(1) > 0`, so the last increment is positive.**
  - Via Lemma FM23 at `r = 3`, the last increment equals `b^2` times a
    positive factor times `R^_8(d,n)`.  `R^_8` coincides with the
    main agent's independent computation of `Psi_3(1)`.
  - With `x = n - (d-5)(d+4)/10`,
    `R^_8 = 8400x^4 + C_2 x^2 + C_1 x + C_0`, where
    - `C_0 = (12/25) d^2(d-1)^2(d+4)^2(d+5)^2`;
    - `C_1 = (16/5) d(d-1)(d+4)(d+5)(2d+3)(2d+5)`;
    - `C_2 = 48(2d^4 + 16d^3 + 22d^2 - 40d - 175)`;
    - `4C_0 C_2 - C_1^2 = (512/25) d^2(d-1)^2(d+4)^2(d+5)^2 (d^4 + 8d^3 - 89d^2 - 420d - 900)`.
  - Cases: `d <= 8` directly (with a finite table for `x < 0`); for
    `d >= 9` the last factor is positive (its largest real root is
    8.815).
- **`R_8(d,n) > 0` for all `d, n >= 0` (PROVED).**  The same method with
  `x = n - (d-4)(d+3)/6`.  The discriminant factor
  `d^4 + 4d^3 - 35d^2 - 78d - 72` has largest root 5.463; the cases
  `d <= 5` are direct.
- **(c) `F(n+2) >= F(0)`: PROVED for `d = 0` (all `n`), for `n <= 1` (all
  `d`, except `(1,0)`), and for `n = 2, 3` (all `d`).**
  - For `n = 2, 3` (main agent): with
    `bracket = C(d+2n,n)(d+2n+1)(d+2n+2) R_8 - A_3 (n+1)(n+2)(n+3) prod_(i=1)^4 (d+n+i)`,
    the bracket is a polynomial in `d` of degree 12 or 13.
  - Its factored forms are `(d+3)(d+4)(d+5)^2(d+6) S_7(d)/2` and
    `(d+3)(d+4)(d+5)(d+6)^2(d+7) S'_7(d)/6`.  Sturm: no root in
    `[0, oo)`, and the value at `d = 0` is positive.
  - Open for `n >= 4`, `d >= 1`; this is covered by the pending
    certificate FM-T1-TR3-CERT at `y = n+2`.
- Checks (main agent): the packet verifier passes (all symbolic
  identities and boundary values), and the discriminant roots were
  confirmed.  The crude bound `C(d+2n,n) >= C(d+2n,4)` leaves 33
  negative coefficients, so it cannot close (c) for `n >= 4`.
- **(c) and (d) for `n = 2, 3`: PROVED exactly (main agent).**
  - For fixed `n`, every `F_(d,n)(x)` is a rational function of `d`.
  - For `n = 2, 3` and every `x <= n+1`:
    `F(n+2) - F(x) = num/den`, and Sturm gives no root of `num` or `den` in
    `[0, oo)`, with a nonnegative value at `d = 0`.
  - So `F(n+2) = max_x F(x)` for all `d` when `n <= 3`, using neptune's
    `n = 0, 1` slices (with the grouped exception `(1,0)`).
  - Verifier: `character_ring_iter/verify_r3_twopower_smalln.py 2 3`.
  - The `r = 3` family is now open only for `n >= 4`.  That is exactly the
    range where the Gaussian-plus-cubic bound has at least 7% slack
    (FM-T1-TR3-CERT pending).
- **(c) and (d) for `4 <= n <= 160` (all `d`) and for `n >= 161`, `d >= n`:
  PROVED (FM-T1-TR3-CERT, luna_max_jupiter; verified by the main agent).**
  - Normalised cubic: `P(w) = 1 + A_1 x + A_2 x^2 + A_3 x^3`, with
    `x = 2w/D`, `D = d+2n+7`, `S = d+n+4`, `Q = R_8`,
    `A_1 = D B_1/(2SQ)`, `A_2 = 2D^2 B_2/(SQ)`, `A_3 = 2D^3 B_3/(SQ)`, and
    explicit `B_1, B_2, B_3`.
  - Main agent: this equals the independent `R(t_0-w)/R(t_0)` in 8,625
    cases.
  - With the binomial bound (ratio `<= e^(-x)`) and `R_8 > 0`, it suffices
    that `P <= e^x` (the case `P < 0` is trivial).
  - Cones: `Q`, `C_1 = 2SQ - D B_1`, `C_2 = SQ - 4D^2 B_2` and
    `C_3 = SQ - 12D^3 B_3` have all coefficients positive after
    `n = 4+v, d = 3n+u` and after `n = 16+v, d = n+u`.  Hence
    `A_1 <= 1`, `A_2 <= 1/2`, `A_3 <= 1/6`, and
    `P <= 1 + x + x^2/2 + x^3/6 <= e^x`.
  - Exact screen for the complementary `d`-ranges (`4 <= n <= 160`):
    `P <= (1+x/K)^K <= e^x` for every `y` in `[1, n+2]`, with
    `K = 512, 128`.  That is 1,406,802 exact Fraction comparisons.
  - Verifier: `character_ring_iter/verify_r3_twopower_cert.py`, re-run
    by the main agent: PASS.
  - **Remaining region of the `r = 3` two-power family:** `n >= 161`,
    `0 <= d < n`.
  - Main-agent data there: `sup_(x>0) P(x) e^(-x) < 1` on a grid
    (`n <= 10^4`, all `d/n`).  `A_1 < 1`; for small `d/n`, `A_2` is about
    `-5.4` to `-5.9` and `A_3` about `1.35` to `1.53`; for
    `d/n >= 0.3` it is coefficientwise.  So `P(x) <= e^x` should hold for
    ALL `x >= 0` (no `y_-` needed).
  - Degree-5 Taylor comparison: dividing by `x`, the requirement
    `min_(x>0) [(1/2 - A_2)/x + x/24 + x^2/120] >= A_3 - 1/6` holds in the
    limit (about 1.53 versus 1.36 at `d = 0`).

### Direction note (main agent, after user review)

The `r = 3` two-power slice was stopped at the user's direction.  The
target is the FULL cone FM3, and per-family certificates (FM54, FM55,
the `r = 3` slice) do not transfer to the next `r` or factor count.  The
partial `r = 3` results above remain valid.  Its only open region was
`n >= 161`, `d < n`, where a degree-5 Taylor comparison held on a
13,013-point grid.  Work now targets SECTOR-level arguments (FM53 is the
model).

### FM-SEC1 (mars): the `r = 1` two-`hat S` sector -- screen and two kills

- **Reformulation (checked).**  Using `hat S_q = 2h_q + 2h_(q-2) - h_(q-1)h_1`,
  the sector is `A_p(h_kappa h_(q-1) h_1) <= 2A_p(h_kappa h_q) + 2A_p(h_kappa h_(q-2))`,
  where `A_p(E) = sum_lambda c_lambda(E) B_p(lambda)`.
  - Mars checked 2,224 cases (`|kappa| <= 10`, `2 <= p, q <= 5`): no
    negative value, and every equality is a parity zero.
  - The main agent independently matched the direct integral in 1,072
    cases.
- **Kill K27 (no one-hat decomposition).**  `hat S_2^2` is not a
  nonnegative combination of FM53-type terms `hat S_r (x) E'` with `E'`
  genuine.
  - The functional `l = (1/4)[U_1U_1] + (1/2)[U_1U_3] - (1/2)[U_2U_2]` on
    `SU(2) x SU(2)` characters is `>= 0` on all 46 admissible generators.
  - But `l(S_2 S_2) = -1`.
- **Kill K28 (no shape-preserving Pieri comparison).**  For `p = q = 2`,
  the Schur coefficient of `hat S_2^perp F_2` at
  `mu_d = (d+2, d+1, d+1, d)` equals `-1` for every `d >= 0`.  Any tableau
  argument must aggregate across shapes.

### Sector programme in the invariant-count (multigraph) model (main agent; user-approved FM3 direction)

Let `I(kappa)` be the number of Sp(4)-invariants in
`(x)_i Sym^(kappa_i) C^4`.  Equivalently, it is the number of loopless
multigraphs on ordered vertices with degree sequence `kappa` and no
3-nesting.  Main agent: the two counts agree in all 219 ordered degree
sequences with at most 5 vertices and total at most 8.

Write `tau_k` for "append a vertex of degree `k`".  Then
`hat S_p = 2 tau_p + 2 tau_(p-2) - tau_(p-1) tau_1`, and the level step is
`(x-y)^2 = 2S_2 + 4 - S_1^2`.  So every FM3 sector is a signed "append
gadgets" inequality for `I`:

- **FM53 (proved, via branching):**
  `I(kappa + {p-1,1}) <= 2 I(kappa + {p}) + 2 I(kappa + {p-2})`.
  - Checked in 695 cases with no nontrivial equality.
  - The factor 2 is needed: coefficient 1 fails in 276 of 695 cases, and
    the maximum of `I(kappa+{p-1,1}) / (I(kappa+{p}) + I(kappa+{p-2}))`
    is 1.81.
- **H-only `r = 2` sector (open):**
  `3 I(kappa + {1,1}) <= 4 I(kappa + {2}) + 8 I(kappa)`, the fifth pass's
  gadget formula.
- **`r = 1`, two `hat S` (open):** the product of two FM53 gadgets.

No Schur shapes appear.  So the shape obstructions K22, K27 and K28 do
not apply to injections in this model.

**Calibration step:** re-prove FM53 by an explicit injection, from these
multigraphs into two copies of each target set.
- **FM-MG1 (saturn): the naive calibration map fails; the inequality stands.**
  - The rule "delete the leaf `v` if `uv` is an edge, otherwise merge `v`
    into `u`" has fibres larger than 2 in 252 of 1,024 cases.  The smallest
    is `kappa = (1,1,1)`, `p = 3`: three graphs merge to one star.
  - The cardinal inequality holds in all 1,024 cases (ordered
    compositions with sum `<= 8`, `p <= 5`).  The maximum ratio is 23/13,
    at `kappa = (1,1,2,2,2)`, `p = 4`.
  - Unused capacity: the second copy of the `p-2` targets.  Next design:
    reroute colliding graphs there by rewiring `{vw, ux} -> wx`.  This
    keeps every degree in `kappa` and lowers `u` to `p-2`.

**Literature check, second pass (main agent, web, 2026-09-27).**
- *"Graphical proof of Ginibre's inequality" (arXiv 2406.08944, J. Stat. Phys. 2024).*
  - A combinatorial proof for XY (U(1)).  It uses the multigraph
    representation of van Engelenburg--Lis (2023): edges coloured red/blue
    and oriented, with divergence constraints `dr = f`, `db = g`.
  - The bijection `{dr = f, db = g} -> {dr = f+g} x {db = -f+g}` reverses
    and recolours the red edges.  It gives an exact product formula, from
    which Ginibre's inequality follows.
  - It depends on abelian currents (orientations and divergences); there is
    no non-abelian extension.  For us it is the template, not a tool: the
    Sp(4) 3-nesting-free multigraph model is the natural non-abelian
    analogue, and a colour/orientation-type involution there is the kind of
    mechanism to look for.
- *Rubey--Westbury (arXiv 1504.02586).*  A combinatorial second
  fundamental theorem for the defining representation of Sp(2n), via
  `(n+1)`-noncrossing perfect matchings.  It covers the defining
  representation only, not products of symmetric powers.  It is the
  noncrossing counterpart of our nesting-free model.
- No SU(2)/Sp(4) Ginibre mechanism was found, consistent with the first
  pass.
- **FM-MG2 (saturn) and the main-agent follow-up: calibration results.**
  - *FM53 (known sector).*  With the vertex order `b_1, u, v, b_2, ...` and
    the local moves (delete `v` if `uv`; merge `v` into `u`, two copies;
    rewire `{vw, ux} -> wx`), a perfect matching exists in all 1,024 cases
    (sum `kappa <= 8`, `p <= 5`).  Main agent: re-verified.
    - Not yet a proof: the matching is algorithmic, and Hall's condition is
      open.
    - The move-graph components reach 524 sources, so Hall's condition is
      not local.
    - Uniform fractional spreading overloads the rewire targets (loads up
      to 8/3).
  - *`r = 2` H-only sector (open) with the analogous local moves* (delete
    `v_1 v_2`; merge `v_1, v_2` into a degree-2 vertex; rewire
    `{v_1 w_1, v_2 w_2} -> w_1 w_2`; plus moves through one extra base
    edge).  Flow test: infeasible in 115 of 255 cases, and still in 59 of
    255 with the richer moves.
  - The failures are near-tight cases of the inequality: slack `L = 2`
    against total `32` at `(1,1,1,3)` and `(2,2,2)`; `L = 10` against
    `280` at `(1,1,1,1,1,3)`.
  - **Conclusion.**  Local injections suffice (computationally) for the
    known sector FM53 but not for the first open sector.  Because the
    open inequality is nearly tight, a proof must be close to an exact
    bijection.  The analogue is the exact product formula in the XY
    graphical proof, not a lossy injection.

### Euler-characteristic mechanism and the graded `r = 2` conjecture qFM3_2 (main agent, 2026-09-27)

Notation: `h_k = chi(Sym^k W)` with `W = C^4`, so `H_(k+1) = h_k`.  Write
`phi_r(kappa) = (1/2) E_(sc x sc)[(x-y)^(2r) prod_i h_(kappa_i)(x,y)]`.
At `r = 2` this is `L(kappa)`.  The H-only sector at level `r` asks
`phi_r(kappa) >= 0`.

**(1) All-ones family: symmetric product (checked, `0 <= r <= 9`, `0 <= m <= 15`).**

    phi_r(1^(2m)) = (2r)! (2r+2)! (2m)! (2m+2)! / ( 4 r!(r+1)! m!(m+1)! (r+m+1)!(r+m+2)! ).

- It is symmetric in `r <-> m`, since `y -> -y` swaps `x+y` and `x-y`.
- At `r = 2`: `L(1^(2m)) = 60 I(1^(2m)) / ((m+3)(m+4))`.
- It is a single family and gives no sector information.

**(2) Positivity on the whole multiplicative H-cone holds only at integer levels.**
- Exact integer arithmetic: `phi_r(kappa) >= 0` for every partition `kappa`
  with `|kappa| <= 16` (915 partitions) and every `1 <= r <= 6`
  (`character_ring_iter/verify_phir_hcone.py 16 1,2,3,4,5,6`).  This includes more than `2r`
  factors, so the statement is stronger than the H-only FM3 sector.
- Gauss quadrature on a 400- and an 800-point semicircle grid: positivity
  fails at every non-integer level tested.
  - At `r = 1/2`, `kappa = (6,4)` gives value/(absolute mass) `= -0.143`.
  - Failures also occur at `r = 0.25, 0.75, 1.5, 2.5`.
  - The integer levels reproduce 0 to `1e-15`.
- So no analytic interpolation in `beta = 2r` can prove the sector.  The
  mechanism is algebraic.
- A Newton (binomial) expansion `phi_r = sum_j a_j C(r,j)` has `a_1 < 0`
  already at `kappa = (2,2)`, so it is not the mechanism either.

**(3) Euler characteristic.**  Let `H = SU(2) x SU(2)`, `P = C^2 (x) C^2`
(so `det(1 - h|P) = (x-y)^2`), and let `M = (x)_i Sym^(kappa_i) W`.
- For ANY action of the abelian Lie algebra `a = P^(+r)` on `M` by
  commuting `H`-equivariant operators, the Koszul complex gives

      2 phi_r(kappa) = sum_k (-1)^k dim Hom_H(Lambda^k a, M).

- So `phi_r >= 0` follows as soon as one such complex has `H`-invariant
  cohomology in even degrees only.
- Slot actions: let `n = Hom(B,A)` sit inside `gl(4)`, where `W = A + B`
  and `A = C^2`, `B = C^2`.  Let `n^(i)` be its action on the `i`-th
  tensor factor, and set `rho(p (x) e_j) = sum_i C_(ji) n^(i)(p)`.  These
  operators commute because `n` is abelian.
- `r = 1`, `C = (1,...,1)`: by Kostant's theorem for the `GL2 x GL2`
  parabolic, the invariant cohomology sits in degrees 0 and 4 only, each
  of dimension `I(kappa)`.  This recovers `phi_1 = I`.
- `r = 2`, generic `C` (Vandermonde; `character_ring_iter/koszul_slot_cohomology.cpp`, rank mod
  `2^31 - 19` over the four weights `(0,0),(2,0),(0,2),(2,2)`):
  - The invariant cohomology is even in all 30 cases finished so far
    (`|kappa| <= 8`).  Example: `(1^6)` gives `[6,0,5,0,18,0,5,0,6]`,
    Euler characteristic `40 = 2 L`.
  - SEE ITEM (17): this fails at size 10 (`(4,1^6)`), so evenness is not
    a general mechanism.
  - The pattern is palindromic in degree.
- `r = 3`, generic `C`: odd cohomology appears.
  - `(1,1)` gives `[1,0,4,0,11,6,8,6,11,0,4,0,1]`.
  - It also appears at `(2,1,1)` and `(1,1,1,1)`.
  - The slot-Koszul route fails at `r = 3` for small `kappa`, and at
    `r = 2` from size 10 (item (17)).

**Relation to earlier Koszul attempts.**
- **K5** tensors `F` with `m` copies of the two-term complex `C^2_x -> C^2_y`,
  using one `n`-action.  It found odd cohomology at `F = S_1^2`, `m = 2, 3`,
  and showed that differentials preserving the `GL(4)`-isotypic splitting
  are forced to fail at odd `m`.
- **K9** shows that differentials natural in the `Sp(4)`-module `M` cannot
  work, so the maps must use the symmetric-power (tensor-factor) structure.
- **FM16** sketched the `r = 2` integral as `chi_y`-genera on `Gr_2(C^4)`.
- The slot action in (3) is new relative to these:
  - It uses `r` copies of `P`, one differential per copy, acting through
    individual tensor factors (as K9 requires).
  - It mixes `GL(4)`-isotypic components, which escapes the K5
    obstruction.
  - Its level is `m = 2r`, which is even.

**(4) Conjecture qFM3_2 (graded `r = 2`).**  Let `U'` be the four
short-root weights `+-e1+-e2` of `Sp(4)`, and let `K_(lambda,kappa)(t)` be
the charge Kostka--Foulkes polynomial.  Then for every partition `kappa`,

    Phi_kappa(t) := sum_(l(lambda) <= 4) K_(lambda,kappa)(t) Ehat_lambda(t)  lies in  N[t],

    Ehat_lambda(t) = < S_lambda(C^4)|_Sp4 , det(1 - t g | U') >_Sp4
                   = [5]_t        if lambda = (m,m,m,m),
                     1 + t^4      if lambda = (a,a,b,b), a > b,
                     t^2          if (lambda_1-lambda_2, lambda_3-lambda_4) in {(2,0),(0,2)},
                     -(t + t^3)   if (lambda_1-lambda_2, lambda_3-lambda_4) = (1,1),
                     0            otherwise.

- The closed form was checked against the direct `Sp(4)` torus integral
  for all 94 shapes with `|lambda| <= 10`.
- At `t = 1`, `Phi_kappa(1) = L(kappa)`, so qFM3_2 implies the whole H-only
  `r = 2` sector.
- **Checked:** all 525 partitions with `|kappa| <= 16` give `Phi_kappa(t)`
  in `N[t]` (`character_ring_iter/verify_qfm3_r2.py 16`, which builds
  `character_ring_iter/kf_charge_sp4.cpp` to enumerate the SSYT with at
  most 4 rows, with charge).
- Equivalent forms:
  - `Phi_kappa(t) = int_Sp4 Q'_kappa(g;t) det(1 - t g|U') dg`, where
    `Q'_kappa` is the modified Hall--Littlewood function.
  - The symmetric function `sum_lambda Ehat_lambda(t) s_lambda` is
    Hall--Littlewood `P(x;t)`-positive.
- **Origin.**  It is the graded Euler characteristic of the `r = 2`
  Koszul complex for `a = n + n t` acting on the FUSION product of the
  `Sym^(kappa_i)`, whose graded multiplicities are Kostka--Foulkes
  polynomials.  Even `H`-invariant cohomology of that complex in each
  degree would prove qFM3_2.  That evenness is FALSE (item (10)); the
  complex is only bigraded-pure (item (16)).
- The pairing is necessary: `sum_lambda K_(lambda,kappa)(1) Ehat_lambda(t)`
  (ungraded Kostka, graded kernel) is NOT positive.  For example, at
  `(2,1,1)` it equals `1 - t + 2t^2 - t^3 + t^4`.
- **No `r = 3` analogue found.**  Replacing `det(1 - tg|U')` by `D_t^2`,
  `D_t D_(t^2)`, `D_t D_1`, `D_t D_(-t)` or `D_t^3` gives negative
  coefficients (for example at `(2,1,1)`).  So do the fusion gradings
  with copy degrees `(0,1,2)`, `(0,1,3)`, `(0,2,3)`, `(0,1,4)`,
  `(0,2,4)`, `(0,3,4)`, `(0,2,5)`, `(0,3,6)`.  At `r = 2`, the degrees
  `(0,2)` and `(0,3)` also fail; only `(0,1)` works.

- **No graded level step.**  Consider graded versions of the level
  recursion,
  `a(t) Phi_(kappa+{2}) + c(t) Phi_kappa - b(t) Phi_(kappa+{1,1})`, with
  `a, b, c` nonnegative Laurent polynomials of degree range `[-6, 6]` and
  `a(1) = 4`, `b(1) = 3`, `c(1) = 8`.  Linear programming (main-agent scratch LP)
  finds no choice that is coefficientwise `>= 0` for all `|kappa| <= 6`.
  So `r = 3` positivity is not a fixed-weight graded consequence of
  qFM3_2.

**(5) FM53 in two lines, via Kostant (main agent).**  Take `r = 1`, one plus
label `p >= 1`, and any polynomial GL(4)-module `E` (not only products of
symmetric powers).  Then

    word = (1/2) int_H S_p chi_E det(1 - h|n)
         = (1/2) sum_lambda c_lambda(E) sum_(w in W^P) (-1)^(l(w)) [ types_w(lambda) in {(p,0),(0,p)} ].

- The second equality uses `det(1-h|n) chi_(V_lambda) = sum_w (-1)^(l(w)) chi_(L(w.lambda))`
  (Kostant / Weyl).
- `S_p = V_p (x) 1 + 1 (x) V_p` pairs only with the types `(p,0)` and
  `(0,p)`.
- The six types are:
  - `l = 0`: `(a,c)`;
  - `l = 1`: `(a+b+1, b+c+1)`;
  - `l = 2`: `(a+b+c+2, b)` and `(b, a+b+c+2)`;
  - `l = 3`: `(b+c+1, a+b+1)`;
  - `l = 4`: `(c,a)`.
- Both odd-length types have both coordinates at least 1, so they never
  match.
- Hence every term is `>= 0`.  This reproves FM53 without FM51 and
  extends it to every GL(4)-module `E`.
- The same bookkeeping explains why two `hat S` factors (or `r >= 2`) are
  harder.  With two factors, `S_p S_q` contains the mixed types `(p,q)` and
  `(q,p)`, which the odd types `(a+b+1, b+c+1)` can hit.  At `r = 2`, the
  generic slot-Koszul cohomology has the types `(0,2)` and `(2,0)` in odd
  degrees (`kappa = (1,1)`).
- The `r = 1` two-`hat S` values are nonnegative in 344 cases
  (`|kappa| <= 8`, `2 <= p <= q <= 5`).  Three tried gradings fail:
  charge Kostka--Foulkes, and that combined with the Kostant length grading
  in either direction (`(-t)^l` or `t^(4-l)`).

**(6) Reduction of `r = 2` evenness to a five-term row complex (main agent).**
Let `a = n + T`, where `T = sum_i c_i n^(i)` (generic slot action, or
the `t`-action on the fusion product).  Use the Hochschild--Serre
spectral sequence for the ideal `n` of `a`:

    E_2^(p,q) = H^p( T ; H^q(n, M) ).

`H^q(n, M)` is given by Kostant.  Its `l(w) = q` pieces have the types
listed in (5).  The `H`-invariant `E_1^(p,q) = Hom_H(Lambda^p P, H^q(n,M))`
matches types against `Lambda^p P`, whose types are:
`(0,0)` for `p = 0, 4`; `(1,1)` for `p = 1, 3`; `(2,0) + (0,2)` for `p = 2`.

- **Rows `q = 1, 2, 3`.**  They meet `Lambda^p P` only for
  `lambda = (m,m,m,m)`, at `(p,q) = (1,1), (3,1), (2,2), (1,3), (3,3)`.
  All of these have even total degree.
- **Row `q = 0` (and its mirror `q = 4`).**  Here `H^0(n,M) = M^n`, and
  the row is

      X_0 -> X_1 -> X_2 -> X_3 -> X_4,
      X_0 = X_4 = (+)_(lambda=(a,a,b,b)) Mult_lambda,
      X_1 = X_3 = (+)_((a,c)=(1,1)) Mult_lambda,
      X_2 = (+)_((a,c) in {(2,0),(0,2)}) Mult_lambda,

  where `dim Mult_lambda = K_(lambda,kappa)`.  The maps come from `T`.
  They move one box from row 3 or 4 to row 1 or 2 (Gaudin-type operators
  on multiplicity spaces).
- **Consequence.**  `H^odd = 0` for the `r = 2` complex follows from
  `H^1(X) = H^3(X) = 0`.
- **Euler characteristics.**  `2 L = 2 chi(X) + 6 sum_(m^4) K`, and
  `chi(X) = 2 sum_(aabb) K - 2 sum_((1,*,1)) K + sum_((2,*,0),(0,*,2)) K`.
  The necessary condition `chi(X) >= 0` holds for all 159 partitions
  with `|kappa| <= 12`.
- **`H^1(X) = 0` as a Poincare lemma for `T`.**  Every `H`-map
  `phi : P -> M^n` with `T(p) phi(p') = T(p') phi(p)` has the form
  `phi = T(.) v`, with `v` in `(M^n)^H`.
- For `kappa = (1^n)`, the natural candidate model of the fusion product is
  `(W^(x)n (x) R_n)^(S_n)`, where `R_n` is the coinvariant algebra, with
  `T = sum_i n^(i) x_i`.  It has the right graded character
  `sum_lambda s_lambda Ktilde_(lambda,1^n)(t)`.  The module isomorphism is
  not yet checked.
- This is the proposed attack on qFM3_2.  The graded version is the same
  row complex, with `t`-graded multiplicity spaces and maps of degree `+1`.
- **Correction (direct computation, main-agent scratch script `rowX.py`).**  The sufficient
  condition `H^1(X) = H^3(X) = 0` is FALSE at `kappa = (1,1,1,1)`, where
  `H(X)_inv = [2,0,1,1,1]`.
  - The full `r = 2` complex there is still even:
    `[2,0,1,0,6,0,1,0,2]`.  So the odd `E_2` class must be cancelled by a
    higher differential from the `(m,m,m,m)` rows `q = 1, 2, 3`.
  - `X` is exact at odd positions for the other 17 partitions with
    `|kappa| <= 6`.
  - A proof along these lines must therefore use the whole spectral
    sequence, not the row `X` alone.
  - `H(X)` is not palindromic, so `X` is not self-dual.

**(7) FM-R3M (luna_max_neptune; main agent re-checked the classification
logic): factor-wise actions at `r = 3`.**
- For `kappa = (1,1,1,1)` at `r = 3`, the `H`-invariant cochain dimensions
  are `[10,48,156,368,690,960,1080,960,690,368,156,48,10]` (Euler
  characteristic 40).  An abstract complex with only even cohomology is
  dimensionally possible, so dimensions give no obstruction.
- **Classification.**  As an `H`-module,
  `gl(4) = (V_00 + V_20) + (V_00 + V_02) + V_11 + V_11`.  So the `H`-maps
  `P -> gl(4)` are spanned by `n_+ = Hom(B,A)` and `n_- = Hom(A,B)`.
  Commutativity on the tensor product must hold factor by factor, and it
  forces each factor to use only `n_+` or only `n_-`.  Hence every action
  through `gl(4)^(+n)` on the tensor factors is an orientation pattern
  `epsilon` in `{+,-}^n` together with a matrix `C`.
- **Random `C`.**  For all-`+` orientation and random `C`, `(1,1,1,1)` at
  `r = 3` gives `[2,0,6,0,10,4,12,4,10,0,6,0,2]`, identical to the
  Vandermonde choice.  At `r = 2`, random and Vandermonde `C` agree on
  `(2,1,1,1,1)`: `[3,0,2,0,6,0,2,0,3]`.
- **Mixed orientations at `r = 2`** (main agent, `kzo.cpp`) do not help.
  - `(1,1)` and `(2,2)` with `(+,-)` give cohomology only in degree 4.
  - `(1,1,1,1)` with `(+,+,-,-)` gives `[1,0,0,1,12,1,0,0,1]` (odd).
  - `(2,1,1)` with `(-,+,+)` gives `[0,0,0,1,6,1,0,0,0]` (odd).
  - All-`+` stays even in every case computed here (sizes `<= 8`).  It
    fails at size 10 (item (17)).
- **Mixed orientations at `r = 3`, `kappa = (1,1,1,1)`** (main agent,
  `kzo.cpp`; random `C`, at least two matrices per pattern, all
  consistent):
  - one `-` factor: `[0,1,1,0,21,3,4,3,21,0,1,1,0]`;
  - two `-` factors: `[1,0,2,7,11,0,26,0,11,7,2,0,1]`;
  - all `+`: `[2,0,6,0,10,4,12,4,10,0,6,0,2]`.
  - With duality (all `-` is the dual of all `+`), every orientation
    pattern has odd invariant cohomology.
  - **Hence no factor-wise action (through `gl(4)^(+4)`) gives even
    cohomology at `r = 3` for `W^(x)4`.**  An `r = 3` Koszul proof would
    need actions that couple tensor factors, or a different module or
    complex.

**(8) The `y -> -y` duality (main agent).**  The semicircle is symmetric, so
`U_n(-y) = (-1)^n U_n(y)` gives, under `y -> -y`:

- `S_p -> S_p` for `p` even, and `S_p -> D_p` for `p` odd;
- `D_q -> D_q` for `q` even, and `D_q -> S_q` for `q` odd.

This maps FM3 words to FM3 words.
- **Example.**  The `r = 2` H-only word with exactly four parts,
  `D_(q_1) D_(q_2) D_(q_3) D_(q_4)` with `q_i = kappa_i + 1`, becomes
  `prod_(q_i odd) S_(q_i) prod_(q_i even) D_(q_i)`.
- **Consequence.**  qFM3_2, restricted to four parts with two odd and two
  even `q_i`, implies the `r = 1` two-`hat S` words
  `(1/2) E[S_p S_q D_(q_1) D_(q_2)]` with `p, q` odd and `q_1, q_2` even.
  For example, `phi_2(2,2,1,1) = (1/2) E[S_3^2 D_2^2] = 4`.
- So the H-only `r = 2` sector already covers part of the open `r = 1`
  two-`hat S` sector.  Even plus labels are not reached this way.

**(9) FM-QR2 (luna_max_saturn; main agent checked): boundary theorems and
an independent screen for qFM3_2.**
- **Proposition 1 (checked).**  Let `A_(lambda,kappa)(q)` be the number of
  SSYT of shape `lambda` and content `kappa` with charge `q`.  Then
  `[t^d] Phi_kappa` is a signed count of tableau-and-shift pairs, with
  `sum_(j=0..4)` on `(m,m,m,m)`, shifts `0` and `4` on `(u,u,v,v)`,
  shift `2` on `(2,*,0)` and `(0,*,2)`, and minus shifts `1` and `3` on
  `(1,*,1)`.
- **Propositions 2--4 (proved).**
  - `Phi_(k) = t^2` for `k = 2`, and `0` for every other `k >= 1`.
  - `Phi_(s,s) = 1 + t^3 + t^4`, `Phi_(s,s-2) = t^2`, and every other
    two-part `Phi` is `0`.
  - `Phi_kappa = 0` whenever `kappa_1 >= sum_(i>=2) kappa_i + 3`.
    Proof: `Ehat != 0` forces `lambda_1 - lambda_2 <= 2` and
    `lambda_1 >= kappa_1`, so `|kappa| >= 2 kappa_1 - 2`.
- **Screen.**  Independent code (Python charge algorithm): all 272
  partitions with `|kappa| <= 12` give coefficients `>= 0`, matching
  `verify_qfm3_r2.py`.
- **Open.**  A uniform injection from the negative pairs to the positive
  pairs.  The simple charge-shifted inequalities
  `A_c(1,*,1) <= B_(c+1) + D_(c+1)` and `A_c(1,*,1) <= C_(c+1)` fail in
  200 and 164 of 294 partitions (main agent).  So the cancellation is
  global, consistent with the higher spectral-sequence differentials in
  (6).

**(10) The fusion Koszul complex is NOT pure (main agent, scratch
`fusion.py`).**
- **Construction.**  `F_kappa = gr M` for the `t`-filtration generated by
  the cyclic vector, at points `c_i = i + 2`.
- **Character check.**  The graded dimensions match the cocharge
  Kostka--Foulkes prediction.  For example, `(2,1,1)` gives
  `[35,45,65,15]` and `(1,1,1,1)` gives `[35,45,65,60,35,15,1]`.
- **Result.**  The invariant cohomology of `a = n + n t`, graded by the
  internal degree `e = (t`-degree of output`) - (#n t` inputs`)`, has odd
  classes for every partition tested with at least three parts.
  Totals over `e`:
  - `(2,1,1)`: `[1,1,2,1,1,0,1,0,1]`;
  - `(1,1,1,1)`: `[3,3,4,2,7,0,1,0,2]`;
  - `(3,2,1)`: `[1,1,2,1,1,0,1,0,1]`;
  - `(2,2,2)`: `[1,2,4,2,1,0,1,0,1]`.
- The one- and two-part cases are even.
- **Why the Euler characteristic is still positive.**  The odd classes
  come in pairs with an even class in the same `e`-block.  For example, at
  `(2,1,1)` the blocks `e = 0` and `e = 2` contribute `1 - 1` each.  So the
  per-`e` Euler characteristics still equal `2 Phi_kappa`.
- **Consequence.**  qFM3_2 is not explained by purity of the fusion
  complex.  The generic ungraded complex is even at these `kappa`,
  consistent with semicontinuity (fusion is a special fibre).  A proof of
  qFM3_2 needs either a different filtration or a cancellation argument.

**(11) FM-R3C (luna_max_neptune): actions coupling two tensor factors at
`r = 3`, `kappa = (1,1,1,1)`.**
- **Space of maps.**  `dim Hom_H(P, End(W (x) W)) = 16`.  The pairwise
  ansatz has a 240-dimensional map space: for each copy of `P`, 8
  one-factor maps plus 12 coupled maps for each of the 6 pairs.
  Commutativity is a quadratic system with 15,246 equivariant components,
  and it was not solved.
- **Conjugation gives nothing new.**  Conjugating a factor-wise action by
  an `H`-equivariant invertible `S` (e.g. `S = I + 2 P_(Lambda^2 A)` on two
  factors) gives genuinely coupled commuting actions.  But they are
  isomorphic complexes, so the cohomology is unchanged (odd).
- **Square-zero lemma (checked by the main agent).**  Suppose `a N` lies
  in `K` and `a K = 0` (all products of action operators vanish).  Rank
  bounds and the self-duality `a_(12-k) = a_k` force
  `dim H_1 + dim H_11 >= c_1 - c_0 = 38`.  So such actions always have odd
  cohomology.
- **Verdict.**  Conditional; the pairwise variety is unresolved.  The main
  agent stops this line: the search space is a large quadratic variety,
  and every structured family found reduces to known cases.

**(12) Sign pattern of the level-`r` Schur weights (main agent).**  Let
`2 w_r(lambda) = sum_w (-1)^(l(w)) g_(r-1)(type_w(lambda))`, where
`g_s(a,c) = sum_j (-1)^j C(2s,j) b(a, 2s-j) b(c, j)` is the multiplicity of
`V_a (x) V_c` in `(x-y)^(2s)`, and `b` is the ballot number.  The types
`(lambda_1-lambda_2, lambda_3-lambda_4)` carrying negative weights are:

| `r` | negative types |
|---|---|
| 2 | `(1,1)` |
| 3 | `(1,1)`, `(1,3)`, `(3,1)` |
| 4 | odd-odd pairs with sum `<= 6`, plus some `(0,0)` shapes (7 types) |
| 6 | 16 types |

The negative region grows with `r`, so no fixed local domination pattern
can serve all levels.

**(13) Rational dependence on `r` (main agent, exact fits).**  For fixed
`kappa` with `|kappa| = 2m`,

    phi_r(kappa) / phi_r(emptyset) = N_kappa(r) / ( prod_(i=2)^(m+1) (r+i) * prod_(i=3)^(m+2) (r+i) ),

with `N_kappa` a polynomial of degree `2m`.  The fit is exact for all 40
partitions with `|kappa| <= 8`, checked at `r = 1..24`.  Here
`phi_r(emptyset) = C_r C_(r+1) / 2`.

- **Examples.**
  - `N_(2) = 2r(r-1)`, `N_(4) = 3r(r-1)^2(r-2)`,
    `N_(6) = 4r(r-1)^2(r-2)^2(r-3)`.
  - `N_(2,1,1) = 24(r^2 - 3r + 12)`, which is positive for all real `r`.
  - `N_(4,1,1) = 12 r(r-1)(3r^2 - 29r + 256)`.
- **Why only integer levels are positive.**  The roots at small integers
  (support-type vanishing) make the numerator negative between them.  This
  explains item (2).
- **No falling-factorial positivity.**  The expansion of `N_kappa` in
  `C(r-1, j)` has negative coefficients for 11 of the 40 partitions (e.g.
  `(4,1,1)`, `(2,2,1,1)`).
- **Schur weights.**  The per-shape weights `w_r(lambda) / phi_r(emptyset)`
  often factor into linear terms, but not always.  For example:
  - `(1,1)`: `-2(r-3)(r+2) / den_1`;
  - `(2,1,1)`: `-(r-1)(r+3)(r+4)(r+6) / den_2`;
  - `(2,2)` contains the irreducible factor `r^2 - 5r + 9`;
  - `(3,3)` contains a cubic factor.

**(14) Twisted-trace form and a failed injection (main agent).**
- **Twisted trace.**  After `y -> -y`,
  `2 phi_r(kappa) = tr(z_M | V) = dim V^+ - dim V^-`.  Here
  `V = (M_kappa (x) W^(x)2r)^H`, `z = (1,-1)` in `H` acts on the `M` part
  only, and `V^+`, `V^-` split `V` by the parity of the `B`-degree of the
  `M` part.  Checked: `dim V^+ - dim V^- = 2 phi_r` in all tests.
- **Target.**  Positivity would follow from a natural `H`-equivariant odd
  operator `V^- -> V^+` that is injective.
- **Tested operator.**  Off-diagonal mixed Casimirs between `M` factors and
  `W` factors, `sum E^(i)_(ab) E^(j)_(ba)` with `a` in `A`, `b` in `B`,
  `i` an `M` factor and `j` a `W` factor.  With uniform weights and with
  generic weights `c_ij` alike, the rank on `V^-` falls well short:

| `kappa`, `r` | `dim V^-` | rank, uniform weights | rank, generic weights |
|---|---|---|---|
| `(2,1,1)`, `r = 2` | 152 | 87 | 126 |
| `(1,1)`, `r = 3` | 280 | — | 240 |

- So no quadratic two-factor operator gives the injection.

**FM-CHK18 (luna_max_mercury, independent checker): items (3), (4), (5),
(8), (9) — all ACCEPTED.**
- **C1 (Euler identity and `r = 1` Kostant).**  Accepted.  Scope as
  stated: the degree-0-and-4-only claim is for the standard diagonal
  action.  An arbitrary action (e.g. the zero action) has invariant
  cochains in other degrees.
- **C2 (closed form of `Ehat`).**  Accepted.  All 53 shapes with
  `|lambda| <= 8` were checked by independent exact constant terms, and
  `Phi_kappa(1) = L(kappa)`.
- **C3 (FM53 via Kostant, extension to all polynomial `GL(4)`-modules).**
  Accepted.  Non-cone tests `S_(2,1)`, `S_(2,1,1)`, `S_(3,1,1)` with
  `p = 1..4` give no negative value.
- **C4 (`y -> -y` duality and the four-part consequence).**  Accepted.
  Independently, `E[S_3^2 D_2^2] = 8`, so `phi_2(2,2,1,1) = 4`.
- **C5 (Propositions 2--4).**  Accepted.

**(15) FM-QR2c (luna_max_saturn): no identification with known positive
families.**
- **Tested families.**  Each was computed on the boundary `kappa` values:
  - type-`C_2` Lusztig q-weight multiplicities after branching;
  - Panyushev's generalized q-multiplicities for the positive short roots
    `{e1-e2, e1+e2}`;
  - Lecouvey's type-C tensor q-multiplicities: these specialize to
    `I(kappa)`, not `L(kappa)`;
  - Shimozono--Weyman parabolic Kostka polynomials for the `2+2`
    parabolic, including the nilpotent-orbit case `X_(2,2)`.
- **Result.**  None equals `Phi_kappa`; most fail already at `kappa = (2)`.
- **Why the known theorems don't apply directly.**  The kernel
  `prod_(alpha in U')(1 - t e^alpha)` is exterior (determinant) type, while
  those theorems use symmetric-algebra factors `prod (1 - t e^alpha)^(-1)`.
  In addition, `U'` contains opposite roots, so it lies in no open
  half-space, violating Panyushev's hypotheses.
- No counterexample to qFM3_2 was found.

**(16) Bigraded purity of the fusion Koszul complex (main agent, scratch
`fusion_s.py`).**
- **Two gradings.**  The fusion complex for `a = n + n t` preserves `e`
  (item (10)) and also `s = (A-degree of output) - (cochain degree)`,
  since `n` and `n t` both raise the `A`-degree by one.
- **Observation.**  In every `(e, s)` block, the invariant cohomology sits
  in a single cohomological degree `k(e, s)`.  Checked for `(2,1,1)`,
  `(1,1,1,1)`, `(3,2,1)`, `(2,2,2)`, `(3,1,1,1)`, `(2,2,1,1)`, with zero
  exceptions.
- **Where the odd classes cancel.**  The odd classes of item (10) are
  whole blocks with odd `k(e, s)`.  They cancel against even blocks with
  the same `e` and a different `s`.
  - Example `(1,1,1,1)`, `e = 2`: `k = 4, 3, 1, 0` at `s = -2, 0, 2, 4`,
    each with dimension 1, so the Euler characteristic in that degree is
    `0`.
  - `k(e, s)` is not a function of `s` alone.
- **qFM3_2 in these terms:**

      for every kappa and e:   sum_s (-1)^(k(e,s)) dim H_(e,s)  >= 0.

  Proving it would need to describe the bigraded cohomology dimensions,
  and a sign-reversing matching of odd blocks with even blocks in the same
  `e`.

**FM-CHK19 (luna_max_mercury, independent re-implementation over `Q`):
items (10) and (13) — ACCEPTED.**
- **D1 (fusion non-purity at `(2,1,1)`).**  Graded dimensions
  `[35,45,65,15]`.  Invariant cohomology by `e`:
  - `e = -2`: `[0,0,0,0,1,0,0,0,1]`;
  - `e = -1`: `[0,0,1,0,0,0,1,0,0]`;
  - `e = 0`: `[0,0,1,1,0,0,0,0,0]`;
  - `e = 2`: `[1,1,0,0,0,0,0,0,0]`.
  - The total is `[1,1,2,1,1,0,1,0,1]`, and the per-`e` Euler
    characteristics match `2(t^4 + t^5)`.
- **D2 (rational fits).**  `N_kappa(r) / den_m(r)` matches the exact values
  at `r = 1..10` for `(2)`, `(2,1,1)`, `(4,1,1)`.  This is a finite check,
  not an all-`r` identity.

**(17) CORRECTION: the generic `r = 2` slot complex is NOT even in
general (main agent, exact character computation, scratch
`ugrade_fast.py`).**
- **Refinement.**  Put `u` on the `A`-degree, i.e. the central `U(1)` of
  `GL(2) x GL(2)`.  Every all-`+` slot action preserves
  `s = (A-degree of output) - k`, so each `s`-block has its own Euler
  characteristic.  These are the coefficients of

      Phi^(r)_kappa(u) = (1/2) E[ prod_i h_(kappa_i)(x,y;u) det(1 - u^(-1) h|P)^r ],
      h_k(x,y;u) = sum_(a+b=k) u^a U_a(x) U_b(y).

- **Check.**  At `u = 1` these reproduce `phi_r` exactly.
- **`r = 2`.**  All coefficients are `>= 0` for `|kappa| <= 8`, but two
  partitions of size 10 are negative (both coefficients shown are `2 phi`):
  - `(4,1,1,1,1,1,1)`:
    `{-8:5, -6:1, -4:10, -2:-1, 4:-1, 6:10, 8:1, 10:5}`;
  - `(3,1,1,1,1,1,1,1)`: negative coefficients at `u^(-2)` and `u^4`.
- **`r = 3`.**  Already `(1,1)` has negative coefficients at `u^(-6)` and
  `u^(-4)`.
- **Consequence.**  EVERY all-`+` slot action (any `C`) has odd invariant
  cohomology at `r = 2` for `kappa = (4,1^6)` and `(3,1^7)`.  The evenness
  observed for `|kappa| <= 8` in item (3) does not persist.  Items
  (3)--(6) remain correct as statements about the cases computed.  The
  "generic `r = 2` evenness" conjecture is withdrawn.
- **Unaffected.**  qFM3_2 is a different (Kostka--Foulkes graded)
  statement, verified to `|kappa| <= 16`.

**(18) FM-QR2d (luna_max_saturn; main agent checked the bookkeeping): the
first page of the fusion spectral sequence.**
- **Setup.**  Filter by the `n t`-cochain degree `q` and take
  `n`-cohomology first (Kostant for the `2+2` parabolic).
- **First page.**  The `E_1` dimensions are
  `a_kappa(p,q;e,s) = sum_lambda Ktilde_(lambda kappa, e+q) * #{I : p(I) = p, Levi type in R_q, mu^I_1 + mu^I_2 - q = s}`.
  Here `I` runs over the six shuffles, and `R_q` is the set of
  `Lambda^q P` types.
- **Euler identity.**  `2 Phi_kappa(t) = t^(N(kappa)) sum_e chi_kappa(e) t^(-e)`,
  with `N(kappa) = sum (i-1) kappa_i`.
- **Match with the data.**  Every reported cohomology block occurs in
  `E_1` with the same dimension.
  - The surplus `E_1` pairs cancel through `d_1`, and at least one `d_2`:
    `E_2^(4,1)(3,-4) -> E_2^(3,3)(3,-4)` at `(1,1,1,1)`.
  - The same `d_2` locations appear at `(2,2,1,1)`.
  - Examples: `(1^4)` goes from 36 dimensions on `E_1` to 22 in `H`;
    `(2,2,1,1)` from 40 to 22.
- **Open.**  Uniform ranks of the higher differentials, or a same-`e`
  matching of odd blocks to even blocks.  Verdict: conditional.

**(19) Cyclage matchings do not give qFM3_2 (main agent, scratch
`cyclage_match.py`).**
- **Setup.**  The negative tableau-shift pairs are `(T, j)` with `T` of
  type `(1,*,1)` and `j` in `{1,3}`.  They should match positive pairs
  with the same `t`-exponent.  Edges are Lascoux--Schutzenberger cyclage
  paths; `charge(cyc T) = charge(T) + 1` was checked in every case.
- **Result.**  Even with up to 8 cyclage steps, and with tableaux of more
  than 4 rows allowed as intermediate nodes, bipartite maximum matching
  leaves 8 of the 40 partitions with `|kappa| <= 8` unmatched.  For example,
  `(1^8)` matches 270 of 280 and `(2,2,2,2)` matches 18 of 20.
- **Consequence.**  The cancellation in qFM3_2 is not local in the cyclage
  graph.

**(20) What `(PCP)` itself needs from the H-only analysis (main agent).**
- By Corollary FM2, `(PCP)` at `m = 2r` minus labels needs `phi_r(kappa) >= 0`
  only for `kappa` with at most `2r` parts `>= 2`.  Parts equal to 1 are
  unlimited, since each plus label `1` gives `S_1 = h_1`.
- The coefficient of `V_a` for `a >= 2` in `G_R` is a one-`hat S` word.
- So `(PCP)` is weaker than the full H-cone statement of item (2).  But it
  still requires every level `r`, and at each level it requires infinite
  families together with the `hat S` words.
- Restricting to `(PCP)` does not remove the `r >= 3` obstruction.

**(21) The graded kernel is a Macdonald weight (main agent).**
- **Reading.**  `det(1 - t g|U') = prod_(alpha short) (1 - t e^alpha)`.
  Together with the `Sp(4)` Weyl density this is the `C_2`
  Macdonald--Koornwinder weight
  `prod_(long) (e^alpha; q)_1 prod_(short) (e^alpha; q)_2` at `q = t`
  (so `k_long = 1`, `k_short = 2`).
- So qFM3_2 pairs the Hall--Littlewood `Q'_kappa(.; t)` with the `C_2`
  Macdonald weight at `q = t`.
- **Level `r`.**  The natural analogue is `(e^alpha; q)_r` on short
  roots, i.e. `prod_(j<r) D_(q^j)`.
- **Test** (`|kappa| <= 10`, 82 partitions): Hall--Littlewood parameter
  `t^b` against weight parameter `q = t^a`.
  - `r = 2`: only `q = t^(+-1)` with `b = 1` is positive.  Every other
    tested coupling fails, in 23--51 cases.
  - `r = 3`: every tested coupling fails.  The least bad, `q = t` with
    `b = 1`, fails in 21 cases (first at `(2,1,1)`).
- So the exact coupling that makes `r = 2` work has no `r = 3` analogue in
  this family.

**(22) FM-R3X-B (luna_max_venus): combinatorial mechanisms for
`r >= 3`.**
- **Fusion-path identity (proved; exact at `r = 2, 3` for `(1,1)`,
  `(2,1,1)`, `(1,1,1,1)`).**

      2 phi_r = sum_(0 <= b_i <= kappa_i) (-1)^(sum b) sum_j C(2r, j)
                mu(kappa - b, 1^j) mu(b, 1^(2r-j)),

  where `mu` counts SU(2) fusion paths.  So `2 phi_r = |C^+| - |C^-|`.
  Example: `(2,1,1)`, `r = 3`: `1520 - 1504`.
- **Survivor.**  A sign-reversing involution on the colored fusion paths.
  This is essentially the original Q3 combinatorics; no map is
  constructed.
- **Kill: Specht-block refinement of the twisted trace.**  Split by
  `S_(2r)`-isotypic blocks of the auxiliary slots; some blocks have
  negative differences.
  - `r = 2`: `d_(2,1,1) = -2` at `(2,1,1)` and `-4` at `(1,1,1,1)`.
  - `r = 3`: `d_(3,2,1) = -2` at `(2,1,1)` and `-4` at `(1,1,1,1)`.
- **Kill: total positivity of the sign-adjusted kernel.**  Take the
  kernel `(-1)^a g_r(a,b)` of `(x-y)^(2r)`.  Its 2x2 minor on labels
  `{0,2}` is `-21` at `r = 2` and `-756` at `r = 3`.

**(23) FM-R3X-A (luna_max_jupiter): representation-theoretic framings for
`r >= 3`.**
- **Two exact reformulations.**
  - Dirac induction for `Spin(5)/Spin(4)`:
    `<1, D-Ind(M_kappa (x) delta^(2r-1))>_G = <M_kappa, delta^(2r)>_H = 2 phi_r`,
    with `delta = chi_(S^+) - chi_(S^-)` and `delta^2 = det(1 - h|P)`.
  - An induced `osp(4|2r)` module, whose odd positive part restricts to
    `P^(+r)`, with invariant superdimension `E - O = 2 phi_r`.
- **Checks.**  `E, O` were computed exactly and agree with item (22):
  - `r = 2`: `(38,32)`, `(156,152)`, `(300,288)`;
  - `r = 3`: `(308,280)`, `(1520,1504)`, `(2792,2752)`.
- **Assessment (main agent).**  Both are restatements of the Euler
  characteristic of item (3), not positivity mechanisms: the needed
  chirality or superdimension inequality is FM3 itself.  Also,
  `osp(4|2r)` is of type II, so the "type-I Kac module" framing does not
  apply literally; the underlying identity with `Lambda(P^(+r))` is
  correct.

**(24) Two cone kills in Proposition 13 coordinates (main agent).**

In the coordinates `z = uv`, `w = u/v` of CENTRAL_CHARACTER_Q3_SEARCH
Proposition 13, an H-only word at level `r` is

    (1/4) sum_(a,b) k_ab(kappa) Delta_r(a,b),
    Delta_r(a,b) = P_2(a) P_0(b) + P_0(a) P_2(b) - 2 P_1(a) P_1(b),
    P_k(a) = CT[C_2(u)^k rho_r(u) u^a],   rho_r = (-1)^r (u - u^(-1))^(2r) >= 0 on |u| = 1.

Here `k_ab >= 0` are the `(u,v)`-coefficients of `prod_i H_(q_i)`.  So `r`
enters only through `Delta_r`, independently of `kappa`.

- **Kill 1: coefficientwise positivity.**  `Delta_r(a,b) < 0` occurs only
  for `a` and `b` both even (e.g. `(2,2)` at `r = 0`), for `r = 0..5`
  with `|a|, |b| <= 12`.  But `k_ab` of H-only kernels is supported
  exactly on that lattice: `a = b = |kappa|` mod 2, and `|kappa|` is even.
  So coefficientwise positivity is unavailable.
- **Kill 2: a log-concavity cone.**  Consider symmetric (including
  `alpha <-> beta`) `SU(2)^2` weight arrays that are log-concave on the
  lattice `{alpha + beta` even`}`.  This cone contains the tents `h_n` and
  excludes `V_1 (x) V_1`.  Of 5,991 random examples, 108 give negative
  values at `r = 1`, and there are negatives at every `r <= 5`.
  - The `r = 1` value is the `Sp(4)`-invariant count of a virtual module,
    so genuine `Sp(4)` (indeed `GL(4)`) structure is essential.  A cone
    defined by weight-shape conditions alone cannot work.

**(25) No positive expansion in standard bases (main agent).**
- **Statement.**  H-only positivity at level `r`, for all `kappa`, says
  that `F_r = sum_(l(lambda)<=4) w_r(lambda) s_lambda` is
  monomial-positive.
- **Sufficient conditions tested.**  A nonnegative expansion of `F_r` in
  `h_mu`, `e_mu` or `p_mu` would suffice.
- **Result.**  None holds.  At `r = 2` and `r = 3`, all three expansions
  have negative coefficients by degree 4, and some already in degree 2.
  - `r = 2`, degree 4: `h`-coefficient `-17` at `(2,1,1)`.
  - `r = 3`, degree 4: `e`-coefficient `-38` at `(3,1)`.
  - The number of negative coefficients grows with the degree.

**(26) A no-go for positivity-only arguments, and the all-ones proof
(main agent).**

Setup: in Proposition 13 coordinates put `c = C_2(u) = 2 cos 2 theta` and
`c' = C_2(v)`, both in `[-2, 2]`.  Then

    (x-y)^2 = (2-c)(2-c'),   xy = c + c',   x^2 + y^2 - 4 = c c',

and the `SU(2)^2` Haar density is `(1/4)(c - c')^2`.  Every H-only kernel
factors as `sum_sigma A_sigma(u) A_sigma(v)`, so

    phi_r(kappa) = sum_sigma (w_sigma / 2) det [[ l_sigma(g), l_sigma(c g) ], [ l_sigma(c g), l_sigma(c^2 g) ]],

where `g = (2-c)^r` and `l_sigma` pushes `A_sigma(e^(i theta)) d theta`
forward to `c`.

- **All-ones family: a short proof, valid for every real `r >= 0`.**  Here
  `h_1 = C_1(u) C_1(v)` is a single feature.  So `l` is the positive
  measure `(2 cos theta)^(2m) |2 sin theta|^(2r) d theta` pushed to `c`,
  and the determinant is `>= 0` by Cauchy--Schwarz.  This matches the
  product formula of item (1).
- **No-go.**  Consider the claim that `int_H chi_M f(c) f(c') >= 0` for
  every polynomial `f >= 0` on `[-2, 2]`.
  - That claim would extend by uniform approximation to
    `f = (2-c)^s`, i.e. to `|x-y|^(2s)`, for non-integer `s`.  That is
    false (item (2)).
  - Hence no argument using only nonnegativity of the level weight (any
    Cauchy--Schwarz or convexity argument with a positive weight) can
    prove the integer levels.
  - A proof must use the polynomial, representation-theoretic nature of
    `(x-y)^(2r) = det(1 - h|P)^r`: the virtual `Sp(4)` character `Q_r`, a
    Koszul complex, and so on.

**(27) Which groups realize level `r` (main agent).**
- **What is needed.**  A compact group `G` with a rank-2 regular
  subtorus on which the roots restrict to `{+-2e_i}` once each and
  `{+-e1+-e2}` exactly `r` times each.  That is `4 + 4r` roots.
- **`r = 1`:** `Sp(4)`.
- **`r = 2`:** `SU(4)`, as in FM16.
- **`r >= 3`:** no simple group works.  At `r = 3` there are 16 roots,
  and no simple group has exactly 16 roots.  Only products such as
  `SU(4) x SO(4)^(r-2)` or `Sp(4) x SO(4)^(r-1)` arise.
- **Consequence.**  The FM16 route (`chi_y`-genera on `Gr_2(C^4)`, with
  Bott-type vanishing) would at `r >= 3` have to be carried out on
  products of flag varieties, summed over a lattice of central
  characters of rank `>= 3`.  That route is open, and is the only one in
  this section that treats every `r` by one geometric mechanism.  Item
  (28) shows its central-character pieces are not individually `>= 0`
  even at `r = 2`, so it would need cancellation across characters.

**(28) The FM16 central-character pieces are not individually
nonnegative (main agent, scratch `grass_tau.py`).**
- **The pieces.**  Expanding the `SU(4)` subtorus integral of FM16 in the
  central character splits `phi_2` into pieces `tau_n`, the Laurent
  coefficients of

      F_kappa(c) = int_H prod_i h_(kappa_i)(cA, c^(-1)B) det(1 - c^2 h|P) det(1 - c^(-2) h|P) dh.

  At `c = 1`, `F_kappa(1) = 2 phi_2(kappa)`.
- **Results.**
  - All coefficients are `>= 0` for the 41 partitions with
    `|kappa| <= 8`.
  - At size 10 they fail for `(4,1^6)`:
    `{... -10:10, -6:-1, 6:-1, 10:10 ...}`, and for `(3,1^7)`.
  - There are 7 failures among the 295 partitions with `|kappa| <= 14`.
- These are the same `kappa`, with the same coefficient pattern, as the
  `A`-degree refinement of item (17).
- **Consequence.**  A per-`n` vanishing argument for the `chi_y`-genera of
  FM16 cannot prove the `r = 2` sector.  The pieces cancel across `n`.
  Of the refinements tested, only qFM3_2's Kostka--Foulkes grading has
  survived every test.

**(29) The bigraded Euler characteristic from characters, and no local
matching (main agent, scratch `bigraded_chi.py`).**
- **Formula.**

      chi_kappa(e,s) = sum_lambda sum_d Ktilde_(lambda,kappa,d) sum_(a,b) (-1)^(a+b) m( Lambda^a P (x) Lambda^b P , V(lambda)_(A-deg = s+a+b) ) [e = d - b]

  Here `m` is the multiplicity of the trivial `H`-type.  It needs only
  Kostka--Foulkes numbers and the `A`-degree split of `V(lambda)`, with no
  linear algebra.
- **Check.**  It reproduces the computed fusion cohomology for `(2,1,1)`
  and `(1,1,1,1)`, block by block, with signs `(-1)^(k(e,s))`.  So under
  bigraded purity it determines the cohomology.
- **Reformulation.**  qFM3_2 is `sum_s chi_kappa(e,s) >= 0` for every `e`.
- **No local matching.**  The negative (odd) blocks pair with positive
  blocks at `Delta s = +2`, `-2` or `-4` depending on the row.
  - Example `(1,1,1,1)`, `e = 2`: the pattern is `+1, -1, -1, +1` at
    `s = -2, 0, 2, 4`, so partial sums in `s` go negative from either
    end.
  - `e = 3`: `+1` at `s = -2`, `-1` at `s = 2`.
  - No monotone or ballot-type rule in `s` proves qFM3_2.

**(30) No safe irreducible functionals beyond `Sym^k` (main agent,
scratch `safe_nu.py`).**
- **Question.**  For which `Sp(4)` irreps `V_nu` is
  `f_nu(M) = <M (x) V_nu, Q_2>` nonnegative on the whole H-cone?
- **Result.**  Tested over `|kappa| <= 7` and `nu_1 <= 5`: only `nu = (k,0)`
  (`Sym^k W`) is.  Every other `nu` goes negative, for example:
  - `f_(1,1)(empty) = -3`;
  - `f_(2,1)((1)) = -2`;
  - `f_(2,2)((1^4)) = -3`;
  - `f_(3,2)((1^7)) = -24`.
- **Consequence.**  An inductive (Pieri) proof of the `r = 2` H-only
  sector cannot close up on a family of individually safe irreducible
  functionals.  An invariant cone would need combinations whose structure
  is unknown.

**(31) No Kostka--Foulkes refinement even for the proven one-`hat S`
sector (main agent).**
- **Test.**  At `r = 1`, write `hat S_p = 2h_p + 2h_(p-2) - h_(p-1) h_1` in
  `GL(4)` terms, and grade each term by `Q'` with shifts:

      2 I_t(kappa+{p}) + 2 t^a I_t(kappa+{p-2}) - t^b I_t(kappa+{p-1,1}),
      I_t(mu) = sum_((a,a,b,b)) K_(lambda,mu)(t).

- **Result.**  No shift pair `(a, b)` in `[-2, 4]^2` is coefficientwise
  `>= 0` over the 40 partitions with `|kappa| <= 8`, for `p = 2, 3, 4`.
  At `t = 1` this is FM53, which is proved.
- **Consequence.**  The graded structure of qFM3_2 is specific to H-only
  words at `r = 2`.  It gives no guidance for the `hat S` sectors.

**(32) Integration-by-parts induction `r -> r+1` is impossible (main
agent, scratch `stein_check.py`).**
- **Stein identity (checked exactly for `r = 0..3`).**  For every
  polynomial `Q`,

      E_r[ (4-x^2)(x-y) Q_x + (2r+1)(4-x^2) Q - 3x(x-y) Q ] = 0,
      E_r = weight (x-y)^(2r) rho(x) rho(y) on [-2,2]^2.

  The same holds with `x` and `y` exchanged.
- **Induction operator.**  So `phi_(r+1) = phi_r o Psi` for
  `Psi F = (x-y)^2 F + T_x[aF] + T_y[a^sigma F]`, where `a` is an arbitrary
  polynomial and `T_x`, `T_y` denote the two null expressions above.
  An induction from `r = 0` would need `Psi` to map products of `h_k` into
  nonnegative combinations of products of `h_k`.  That would also explain
  integer-only positivity, since the base case fails at non-integer `r`.
- **Obstruction.**  Every such null expression vanishes at
  `(x,y) = (2,2)`: each term carries `(4 - x^2)` or `(x - y)`.  So `Psi F`
  vanishes at `(2,2)` for every `F`.  But every nonzero element of the
  H-cone is strictly positive there, since `h_mu(2,2)` is a dimension.
- **Enlarging the cone does not help.**  Adding sums of squares fails,
  because the Leibniz rule would need `phi_r(h_kappa G^2) >= 0`.  That is
  false, since `h_kappa` changes sign.
- **Consequence.**  No local integration-by-parts induction on the level
  can prove the H-only sector.

**(33) Tail positivity of the central-character refinement: true for
`r <= 5`, killed from `r = 6` (main agent, scratch `grass_tau.py`,
`tails2.py`, `tailmap.py`, `wgrade.py`, `adams_check.py`).**
- **Setup.**  For exponents `e_1, ..., e_r` in `{+2, -2}` put

      F^(r)_kappa(c) = int_H prod_i h_(kappa_i)(cA, c^(-1)B) prod_j det(1 - c^(e_j) h|P) dh = sum_n tau_n c^n.

  Then `F^(r)_kappa(1) = 2 phi_r(kappa)`.  Item (28) is the case `r = 2`,
  `(e_1, e_2) = (2, -2)`.
- **The sign pattern only shifts.**  `P` is self-dual and `det(h|P) = 1`,
  so `det(1 - c^(-2) h|P) = c^(-8) det(1 - c^2 h|P)`.  There is one
  refinement per level.  After centring it is palindromic (swap `A`, `B`),
  and only powers `c^n` with `n = N mod 4` occur.
- **Meaning.**  In `Spin(6) = SU(4)` with Levi `L = S(U(2) x U(2))`, the
  spaces `c^(+-2) P` are the two nilradicals of the parabolic of the Klein
  quadric `Q^4 = Gr(2,4)`, and `c^2 = e^(e_1)`.  Other central weights fail
  at once: `(1,-1)` and `(4,-4)` at `|kappa| = 2`, `(2,-2,0)` at size 8.
- **`r = 2` identity.**  Weyl integration relative to `L` gives

      tau_(2m)(kappa) = < h_kappa[W], psi^m(Lambda^2 W) >_(SU(4)),     psi^0 := 6,

  where `psi^m` is the Adams operation and the `W`-orbit of `e_1` has 6
  elements.  This was checked exactly for 12 partitions at every `m`
  (`adams_check.py`, 0 mismatches).  Hence

      2 phi_2(kappa) = sum_(m in Z) < h_kappa[W], psi^m(Lambda^2 W) >,

  which is the `h_kappa`-weighted eigenvalue density of Haar `SO(6)` at
  eigenvalue `1`.
- **General level (derived by the same argument, not separately
  checked).**

      tau^(r)_(2m) = < M, sum_(w in W/W_L) w( e^(m e_1) prod_(alpha in n) (1 - e^(-alpha))^(r-2) ) >_(SU(4)).

  By Atiyah--Bott this is `< M, chi(Q^4, O(m) (x) lambda_(-1)(Omega^(+(r-1)))) >`.
- **Conjecture `T_r` (tail positivity).**  For every `kappa` and every
  `n_0`, `sum_(n >= n_0) tau_n >= 0`.  By palindromy the lower tails
  follow.
  - The case `n_0 = min` is H-cone positivity at level `r`, so `T_r`
    implies FM3's H-only sector at level `r`.
  - Equivalent forms: `F(c) / (1 - c^(-4))` has nonnegative coefficients
    in `c^(-1)`; or there is an injection from negative units to positive
    units of weakly larger `c`-weight.
- **Evidence (no failure anywhere).**
  - `r = 2`: 526 partitions, `|kappa| <= 16`.
  - `r = 3`: 295 partitions, `|kappa| <= 14`.
  - `r = 4`: 160 partitions, `|kappa| <= 12`.
  - `r = 5`: 83 partitions, `|kappa| <= 10`.
  - Also `r = 1` (both signs) to size 14.
  - The individual `tau_n` do go negative (item (28)).
- **How strong it is.**
  - At `r = 2` tails are nearly automatic.  Only 18 of 526 partitions
    have any dip, and the deepest keeps 87.5% of the running maximum, at
    `(4,1^6)`.
  - At `r = 3` and `r = 4` the dips are deep.  At `(2,1,1)` the tail
    falls to 33% (`r = 3`) and 20% (`r = 4`) of its running maximum.
  - So for `r >= 3` it is a genuine strengthening.  The negative mass sits
    at middle `c`-weights and is paid for from the extremes.
- **KILLED for `r >= 6`.**  `T_r` holds at `r = 2..5`: sizes 22, 18, 16
  and 14, that is 2,540, 911, 526 and 295 partitions, with no failure.
  But it fails at `r = 6` from size 6 on.
  - Example `r = 6`, `kappa = (2,2,1,1)`: the upper tails from the top are
    `2, 63, 643, 2668, 5274, 5864, 4750, 2860, 970, -144, 446, ...`, with
    total `5720`.
  - Failures among the 41 partitions with `|kappa| <= 8`: 3 at `r = 6`, 6
    at `r = 7`, 13 at `r = 8`, 15 at `r = 9`, 17 at `r = 10`.  The first
    at `r = 7` is `(2,1,1)`.
  - As `r` grows the Koszul degree dominates the `c`-weight, and partial
    Koszul Euler characteristics alternate.
- **The grading is forced but not enough.**  Take a general grading
  `a (M-weight) + e (per P-factor)` and test `r = 2..9`, `|kappa| <= 8`.
  Only `(a, e) = (1, 2)`, the geometric one, survives at any level, and
  only through `r = 5`.  The choices `(1,0)`, `(0,2)`, `(2,2)`, `(1,1)`,
  `(3,2)`, `(1,4)`, `(2,1)` and `(4,1)` fail at `r = 2` on `kappa = (2)`.
- **Consequence.**  Tail positivity is a low-level phenomenon (`r <= 5`),
  not a uniform mechanism.  It gives nothing for the full cone beyond the
  levels where H-cone positivity is already checked directly.

**(34) Real-level structure: the level polynomial and Conjecture LP (main
agent, scratch `onepart.py`, `levelpoly*.py`, `lpscan.py`, `xcheck.py`,
`wallach.py`, `genG.py`, `adamsG.py`).**
- **Non-integer levels fail at every level, through one-part words.**  The
  exact closed form, fitted to 10 digits for `m = 1..4` at seven real `r`
  and matching every integer value, is

      phi_r((2m)) = (m+1) Gamma(r) Gamma(2r+1) Gamma(2r+3) / [ 2 Gamma(r+2) Gamma(r+m+2) Gamma(r+m+3) Gamma(r-m) Gamma(r-m+1) ].

  - So `phi_r((2m)) = 0` at the integers `r <= m`, and it is `< 0` on
    every open interval `(j, j+1)` with `j < m`.
  - Gauss quadrature counts of negatives among `|kappa| <= 12` are 45, 26,
    14, 7, 3, 1, 0 at `r = 0.5, 1.5, ..., 6.5`.  The negatives all have a
    large first part.
- **Rational reduction (exact).**
  - Put `x = 2cos a`, `y = 2cos b`, `u = (a+b)/2`, `v = (a-b)/2`.  The
    weight is `|2 sin u|^(2r) |2 sin v|^(2r)`, and the Fourier mode
    `e^(2iju)` has normalized mean

        rho_j(r) = (-1)^j r(r-1)...(r-j+1) / ((r+1)...(r+j)).

  - Hence, for every FM3 word `w` (`h_kappa` times `hat S_p` factors),
    `phi_r(w) / phi_r(empty) = P_w(r) / Q_w(r)`, where `P_w` is a
    polynomial and `Q_w = prod (r+i)^(e_i) > 0` for `r >= 0`.
  - Checked against exact Catalan-moment evaluation on 395 (word, level)
    pairs, with 0 mismatches.
  - **FM3 at all levels for `w` is equivalent to `P_w(n) >= 0` for every
    integer `n >= 1`.**
- **Observed shape.**

      P_w(r) = c_w r^e (r-1)^2 (r-2)^2 ... (r-d+1)^2 (r-d) G_w(r),     c_w > 0.

  - The positive integer roots are exactly `1..d`: double, except `d`,
    which is simple.  This is the support-region vanishing
    (`d = kappa_1 - |kappa|/2` for H-only words), with a double zero.
  - Examples, as `phi_r / phi_r(empty)`:
    - `(2)`: `2r(r-1)/((r+2)(r+3))`;
    - `(3,1)`: `60r(r-1)/((r+2)(r+3)^2(r+4))`;
    - `(4)`: `3r(r-1)^2(r-2)/((r+2)(r+3)^2(r+4))`;
    - `(2,2)`: `4(r^4+5r^2+54)/((r+2)(r+3)^2(r+4))`;
    - `(2,1,1)`: `24(r^2-3r+12)/((r+2)(r+3)^2(r+4))`.
- **Conjecture LP.**  For every FM3 word `w`, `G_w(r) > 0` for all real
  `r >= 1`.
  - LP implies FM3 for `w` at every level at once: the value is zero for
    `r <= d` and positive after.  This includes the `hat S` sectors.
  - The no-go of item (26) does not apply, because LP is a statement about
    real `r`.  After the forced integer zeros are divided out, the
    remainder is positive at non-integer levels as well.
- **Evidence.**
  - H-only, `|kappa| <= 10` (83 partitions): `G_kappa` has no positive
    real root at all.
  - Words with `hat S` factors, total degree `<= 10` (167 words, 12 plus
    patterns): the same, with the canonical integer-root pattern in every
    case.
  - First stray root: `kappa = (6,1^6)`, size 12.  `G` has degree 6 with
    real roots `0.537` and `0.766`, so it dips inside `(0,1)`, but it is
    positive on `[1, oo)`.  Hence LP is stated on `[1, oo)`.
  - Extended scan (`lpscan.py 12`): 407 words, namely H-only
    `|kappa| <= 12` and plus labels `(2)`, `(3)`, `(4)`, `(2,2)`, `(3,2)`
    to total degree 12.  No violation, and the only stray root is
    `(6,1^6)`.
  - H-only size 14 (`lp14.log`): 4 stray pairs, all inside `(0,1)`:
    `(7,2,2,1^3)`, `(7,2,1^5)`, `(7,1^7)` and `(6,1^6)`.
  - Boundary partitions `kappa_1 = |kappa|/2` of size 18 (all 30) and
    hooks: stray roots again only inside `(0,1)`.
- **Independent check (FM-LP-X, luna_max_mars).**
  - The rational reduction was re-derived through the Beta recursion
    `(r+j+1) I_(j+1) + (r-j) I_j = 0`, with 0 mismatches against exact
    moments on 50 random (word, level) pairs.
  - An exact Sturm scan covered 844 words: all H-only partitions of sizes
    14 and 16, and `hat S` words with plus labels `(2)`, `(3)`, `(2,2)`,
    `(4,2)`, `(3,3)`, `(2,2,2)` to total degree 14.  In every case the
    zero pattern is canonical, `c_w > 0`, and `G_w` has no root in
    `[0.99, oo)`.
  - It also found the identity `hat S_p(x,y) = U_p(x) + U_p(y)`.
- **KILLED: Conjecture LP is false (hooks, from `k = 14`).**
  - `kappa = (14,1^14)`: `G` has positive real roots `0.00978`, `0.99961`,
    `1.60500` and `1.78378`.  So `phi_r < 0` on `(1.605, 1.784)`, with
    `phi_(1.7)/phi_(1.7)(empty) = -2.60`.
  - The integer values are positive: `1`, `13.6`, `16` at `r = 1, 2, 3`.
  - `(13,1^13)` is still positive at `r = 1.7` (`+13.5`).
  - Found by FM-LP-B (luna_max_jupiter) and FM-LP-X2 (luna_max_mars)
    independently; re-checked by the main agent.
  - FM-LP-B also killed four certificate forms for `G_w`, all at degree
    `<= 4`:
    - coefficient positivity of `G_w(1+s)` in the monomial,
      rising-factorial and binomial bases (first failure `(2,2)` or
      `(2,1,1)`);
    - nonnegative combinations of shifted Fourier atoms
      `rho_i(r-1) rho_j(r-1)`, refuted at `(2)` by an exact separating
      functional;
    - a positive-semidefinite Fourier Gram matrix, which fails for the
      empty word (`v^T C v = -8`);
    - pointwise positivity of the word itself (`h_2(0, 1/2) = -7/4`).
  - Two-part words have an exact finite alternating sum through item
    (35).  `G` is even with nonnegative coefficients on all 28 pairs with
    `alpha + beta <= 14`.  No uniform proof is known.  FM-LP-X2 reports the
    same pair of roots in `(1,2)` for all hooks `k = 14..20`.
  - By contrast, all 72 H-only boundary words of sizes 18 and 20, and 63
    boundary words with one `hat S` factor, have no root in `[1, oo)`.
  - Exact hook values: `phi_1((k,1^k)) = 1` and
    `phi_2((k,1^k)) = (k^2 - 5k + 10)/2`.
  - Exact check: `G_(14,1^14)(17/10) < 0`.
  - The near-1 stray root stays below 1 through `k = 20`.  A fit gives
    `G(1)/G'(1) ~ 2^(-k) (k/4 + 2)`.
  - **Consequence.**  Positivity of `phi_r` is an integer-level phenomenon
    even above `r = 1`.  No real-parameter argument (LP, Wallach-type
    continuous part) can prove the cone.  The exact statements at integer
    `r` survive: the rational reduction, the closed form of item (35), the
    forced zeros, the stability of item (36), and the local transport.
- **Stray roots occur only for boundary words `kappa_1 = |kappa|/2` (so
  `d = 0`).**
  - Hooks `(k,1^k)`: the upper stray root is 0.766, 0.932, 0.972, 0.987,
    0.994 for `k = 6..10`.  It tends to `1` from below, where
    `phi_1 = 1`.
  - Boundary words behave as degenerations of the `d = 1` pattern
    `r(r-1)`.  This is where LP is tight.
- **Two-part `kappa`.**  `G_kappa` is an even polynomial with positive
  coefficients, e.g. `(4,4)`:
  `9r^8 + 426r^6 + 21441r^4 + 64524r^2 + 216000`.  So LP holds manifestly
  there (checked `|kappa| <= 8`).  In general `G` is not
  coefficient-positive, e.g. `(2,1,1)`: `24(r^2 - 3r + 12)`.
- **Generalized weights (kills and one survivor).**
  - The weight `prod_i |det(1 - zeta_i h|P)|^2` with unimodular
    `zeta_i` not `+-1` fails.  For example `Re zeta = 4/5` fails at
    `(4,1,1)`.
  - Adams twists `(psi^m x - psi^m y)^2` fail for `m >= 3` at
    `kappa = (2)`.
  - The mixed weights `(x-y)^(2a) (x+y)^(2b)` survive: 7 pairs `(a,b)`,
    160 partitions, no failure.
- **What this gives the full cone.**  One real-parameter positivity
  statement per word, LP, implies every level of FM3 at once.  Two things
  remain open:
  - the vanishing pattern (support region with double zeros), probably
    from the Fourier support of `w` in `(u, v)`;
  - positivity of `G_w` on `[1, oo)`, for example through a
    positive-weight representation of `G_w`.

**(35) Closed form of the FM12 weights, a real-level product formula for
irreducible characters, and the forced zeros (main agent, scratch
`sp4irr.py`; walk check inline).**
- **Closed form (checked).**  Put `p = (a+b)/2` and `q = (a-b)/2` for an
  `Sp(4)` weight `mu = (a,b)`, `a+b` even.  Then

      N_(2r-1)(a+1, b) = (a+2)(b+1) C(2r+1, r-q) C(2r+1, r-1-p) / (2r(2r+1)).

  - This matches the quarter-plane walk counts exactly for all 240 cases
    with `r = 1..8`, `a < 2r`.  It is a Guy--Krattenthaler--Sagan type
    formula (reflection/LGV).
  - With Lemma FM12,

        Q_r = sum_mu (-1)^b (a+2)(b+1)/(2r(2r+1)) C(2r+1, r-q) C(2r+1, r-1-p) V_mu.

    This reproduces the displayed `r = 2, 3, 4` weight vectors (for
    example `294, -378, 168, 189, -105, -35, 27, 21, -7, 1` at `r = 4`).
- **Real-level product formula (proved from the above).**

      phi_r(V_(a,b)) / phi_r(empty) = (-1)^b (a+2)(b+1)/2 * r(r-1)...(r-q+1) * (r-1)(r-2)...(r-p) / [ (r+2)...(r+q+1) * (r+3)...(r+p+2) ].

  - Both sides are rational functions of `r` (item (34)) that agree at
    every integer `r >= 1`, so they are equal.
  - Checked directly for 15 irreducibles `(a,b)` with `a+b <= 8`.
  - The one-part formula of item (34) is the case `(2m, 0)`.
- **The forced zeros, explained.**
  - For every `mu` in a word, `phi_r(V_mu)` carries
    `r(r-1)...(r-q+1) (r-1)...(r-p)`.
  - For an H-only word with `d = kappa_1 - |kappa|/2 >= 1`, every
    `Sp(4)`-constituent `mu` of `h_kappa` has `p >= d` and `q >= d`.
    Brauer--Klimyk: `mu = (kappa_1, 0) + nu'` with `|nu'_1| + |nu'_2|` at
    most the rest of the degree, and reflections keep the bounds.
    Checked directly (`pqcheck.py`, Racah--Brauer peeling) for all 159
    `kappa` with `|kappa| <= 12`, with no exception.
  - Hence `r(r-1)^2 ... (r-d+1)^2 (r-d)` divides `P_w`.  This is the
    observed pattern of item (34), now derived.  The Brauer--Klimyk bound
    with reflections still has to be written out.
- **Full-cone screen through the closed form (main agent,
  `fm3closed.py`).**
  - The screen computes `m_mu(w)` by Racah--Brauer peeling of the weight
    multiset (with `hat S_p` weights), then sums the closed form.
  - Results: 14,251 FM3 words (H-only and 22 plus-label patterns up to
    five `hat S` factors), total degree `<= 22`, levels `r = 1..20`.
    That is 285,020 (word, level) pairs, with **0 negative values**.
  - Validation (`xcheck_hi.py`): `phi_r(w) = (closed form)/(2r(2r+1))`
    agrees exactly with direct Catalan-moment integration on 28 checks at
    `r = 9, 12, 15, 20`, covering H-only words and words with one or two
    `hat S` factors.
  - This extends the earlier `hat S` evidence (item MP: `r <= 3`, size
    `<= 11`) by far.
- **FM3 in closed form.**  At level `r`, FM3 for a word `w` with
  `Sp(4)` multiplicities `m_mu(w)` reads

      sum_mu (-1)^b (a+2)(b+1) C(2r+1, r-q) C(2r+1, r-1-p) m_mu(w) >= 0.

  In `i = r-q`, `j = r-1-p` (`0 <= j < i <= r`) the weight is
  `(-1)^(i-j-1) (i-j)(2r+1-i-j) C(2r+1,i) C(2r+1,j)`.  This has the shape
  of a non-intersecting pair of lattice paths (LGV), which is a natural
  target for a sign-reversing involution.
- **Local transport (observed).**  In `(p,q)` coordinates the diagonal
  neighbours `(a+-1, b+-1)` are the four grid neighbours, and the sign is
  `(-1)^(p+q)`.  So FM3 at level `r` says the checkerboard function
  `F(p,q) = (-1)^(p+q) W_r(p,q) m_(p,q)(w)` has nonnegative sum.
  - For all 159 `kappa` with `|kappa| <= 12` and `r = 1..7` (1,113 cases),
    the odd-cell mass can be transported to adjacent even cells without
    exceeding their mass.  This is linear-programming feasibility, and
    Hall's condition holds with grid-neighbour moves only.
  - Single-neighbour and two-neighbour greedy bounds fail (for example
    `(1,1)` at `r >= 2`).  So an explicit transport rule has to share each
    even cell among up to four odd neighbours.
  - An explicit rule would be an injection proof of the H-only sector at
    every level.
  - Extended (`transport2.py`): H-only `|kappa| <= 12`, `r <= 9` (1,440
    pairs) all feasible.  Words with plus labels `(2)` and `(2,2)` are
    also feasible with grid-neighbour moves.
  - Words with `hat S_3`, `hat S_4`, `hat S_5` are not.  For example
    `h_1 hat S_3` at `r >= 3` has negative mass at `(a,b) = (2,2)` with no
    positive neighbour.  Allowing moves of size up to 4 in `(a,b)` repairs
    `hat S_3`, `hat S_4` and `hat S_3 hat S_2`, but not `hat S_5`: the
    radius grows with the label.
  - So the `hat S` sectors need a transport organized by
    `hat S_p = V_(p,0) - V_(p-1,1) + V_(p-2,0)` itself, as FM53 does at
    `r = 1`.
  - **Exact confirmation (`exactscan.py`, integer max-flow).**  At high
    `r` the floating-point linear program is unreliable in both
    directions.  One float "infeasible" case, `(4,2,2,1^8)` at `r = 12`,
    is exactly feasible.
    - Exact integer max-flow: all 6,456 (word, level) pairs are feasible.
      These are all H-only `|kappa| <= 16` (even) at `r <= 12`, plus hooks
      `(k,1^k)` for `k <= 20`.
    - The margin is tiny.  The ratio (FM3 total)/(negative mass) drops to
      `1.5e-9`, at `(2,2,2,1^10)`, `r = 12`.  At each size up to 12 the
      tightest word is the all-ones family, which item (26) proves.
  - **Explicit greedy rule (FM-TR, luna_max_saturn).**  Take the odd
    cells by increasing `p` (ties: decreasing `q`).  Each sends its demand
    to `(p-1,q)`, `(p,q+1)`, `(p,q-1)`, `(p+1,q)` in that order, taking
    whatever capacity remains.
    - Independent main-agent check (`greedy.py`): it succeeds for all
      30,740 cases with H-only `|kappa| <= 20` (even) and `r <= 20`.  The
      minimum residual margin is 40.
    - KILLED on hooks: it fails from `(13,1^13)`, `r = 4` (margin
      `-144,144`), with 103 failures among hooks `k <= 24`.
    - In every failing case exact max-flow is still feasible
      (`hookflow.py`: hooks `k <= 24`, `r <= 20`), and FM3 is positive.
    - Search over 96 greedy rules (8 source orders x 24 neighbour orders;
      `greedysearch.py`): none works.  The best, sources by `a = p+q` with
      order `(p-1, q-1, q+1, p+1)`, fails 6 times over 538 words x 16
      levels.
    - So the transport exists in all exact tests but needs a global
      matching, not a one-pass rule.
    - Saturn's Hall-set screen: minimum slack is positive (40) in all
      requested cases, and a minimal-slack set is always an interval in
      source order (staircase bands).
  - All four grid directions are needed (`dirflow.py`).  Every
    restriction to two or three directions in `(p,q)` fails; the best,
    `{p-1, p+1, q+1}`, fails in 43 of 1,440 cases, first at
    `(4,1^4)`, `r = 7`.  So there is no one-dimensional ballot rule.
    Confirmed by exact integer max-flow (`exactdir.py`), with the same
    counts and 0 failures for all four directions.
  - **FM-LP-I (luna_max_venus).**  The concatenated walk model (Pieri path
    `0 -> mu`, then a unit walk `(mu_1+1, mu_2) -> 0` of length `2r-1`,
    sign `(-1)^(mu_2)`) matches exact moments on all 67 partitions with
    `|kappa| <= 8` at `r = 2, 3`.  KILLED: any involution that changes
    only the last Pieri edge and the last tail step.  Fibers with more
    negatives than positives occur in 22 of 67 words at `r = 2` (first
    `(2,1,1)`) and 26 of 67 at `r = 3` (first `(1,1)`).  A working move
    must reach earlier Pieri history, which is consistent with the
    four-direction transport above.
- **Literature (FM-LP-L, luna_max_neptune).**
  - The real-level product formula also follows from Kadell's Theorem 1
    ("An integral for the product of two Selberg-Jack symmetric
    polynomials", Compositio Math. 87 (1993)), the Selberg integral with
    one Schur insertion.  Put `xi, eta = (1 - cos(theta +- phi))/2`; the
    weight is `z^(r-1/2) (1-z)^(-1/2)` with Vandermonde exponent 1.
    Checks: `E[xi+eta] = (2r+3)/(r+3)` and
    `E[xi eta] = (2r+3)(2r+1)/(4(r+3)(r+2))`.
  - The packet matched the weight to Remling--Rosler's compact real
    Grassmannians.  The main agent has not verified this, and it looks off
    by a convention.  Those spaces have long-root multiplicity `1/2`, or
    an extra `2e_i` root in the complex case, whereas here
    `(k_long, k_short) = (1, r)` exactly.
  - No known theorem gives positivity on the word cone.  Verdict: no
    progress on the positivity half.
- **KILL: Heckman--Opdam expansion positivity (hypergroup route; main
  agent, `hoexp.py`).**  Expand words in the Jacobi / Heckman--Opdam
  polynomials `P_lambda` of `Sp(4)` at `k = (1, r)`, built by
  Gram--Schmidt on Weyl-orbit sums with exact moments.
  - At `r = 1` these are the Weyl characters, and all 66 H-only words with
    `|kappa| <= 8` have nonnegative coefficients.
  - At `r = 2` already the generator `h_2` has coefficient `-1/6` on
    `P_(1,1)`.  Negative coefficients occur in 23 of 66 words at `r = 2`
    and 36 of 66 at `r = 3`.
  - So no positive-linearization argument (generators expand positively,
    products stay positive) can prove FM3 at `r >= 2`.
- **PROVED (FM-P35, luna_max_venus; reviewed by the main agent).**
  - **(A) Closed form.**  Guy--Krattenthaler--Sagan, "Lattice paths,
    reflections, & dimension-changing bijections", Ars Combin. 34 (1992)
    3--15, give the quadrant reflection formula

        N_n(A,B) = C(n,K)C(n,L) - C(n,K+1)C(n,L-1) - C(n,K+1)C(n,L+1) + C(n,K+2)C(n,L),
        K = (n+A+B)/2,  L = (n+A-B)/2.

    With `n = 2r-1` and `(A,B) = (a+1, b)`, algebra gives the closed form.
    Exact checks: 240 endpoints, compared with dynamic programming.
  - **(B) Real-level product formula.**  Put
    `phi_r(empty) = Gamma(2r+1) Gamma(2r+3) / (2 Gamma(r+1) Gamma(r+2)^2 Gamma(r+3))`.
    FM12 and (A) give the formula at every integer `r`.  Both sides are
    rational in `r` (item (34)), so they agree for all real `r`.
    Exact checks: 620 cases.
  - **(C) Forced zeros for every FM3 word.**
    - Define `d(w)` as in item (36), the minimum over the
      `hat S_p = 2h_p + 2h_(p-2) - h_(p-1) h_1` expansions.
    - In `(P,Q) = ((a+b)/2, (a-b)/2)` coordinates the `Sp(4)` Weyl group
      acts by signed permutations and `rho = (3/2, 1/2)`.
    - Every weight of the remaining factors has `|P|, |Q| <= N/2`.
      Brauer--Klimyk then gives `p(mu), q(mu) >= d(w)` for every
      constituent.
    - Hence `r(r-1)^2...(r-d+1)^2(r-d)` divides `P_w`, the denominators
      being products of `r+i` with `i > 0`.
    - Exact scan: 299 words of total degree `<= 10`, 0 mismatches.
- **What this gives the full cone.**
  - LP (item (34)) becomes an explicit statement: positivity on `[1, oo)`
    of `sum_mu (-1)^b m_mu(w) c_mu` times products of linear factors.
  - The forced-zero half of LP is reduced to a Brauer--Klimyk bound.
  - The positivity half is the open core.

**(36) Forced zeros for `hat S` words, and reduction to boundary words
(FM-LP-A, luna_max_saturn; main-agent checks `reduce_d.py`, stability
check).**
- **Defect for `hat S` words.**  Expand each
  `hat S_p = 2h_p + 2h_(p-2) - h_(p-1) h_1` into label multisets `gamma`,
  and put `delta(gamma) = max(0, gamma_1 - |gamma|/2)`.  Define
  `d(w) = min` of `delta` over all choices.
  - Example: `hat S_4` gives `d = 1`, matching its simple root at `r = 1`.
  - `d(w)` reproduces the zero pattern of all 299 words of total degree
    `<= 10`.
- **Value zeros below `d` (proved).**
  - `SU(2)` triangle inequalities give the Fourier support
    `C_w(j,j') = 0` whenever `max(|j|,|j'|) < d(w)`.
  - Since `rho_j(k) = 0` for `|j| > k`, this gives `P_w(k) = 0` for
    integers `0 <= k < d(w)`.
- **Double zeros at `1..d-1` and the zero at `d`.**
  - On the Fourier side these are two explicit cancellation identities
    (derivative and endpoint), unproved there.
  - Item (35) gives them directly: each `Sp(4)` constituent contributes
    `r(r-1)...(r-q+1) (r-1)...(r-p)` with `p, q >= d`.  The branching
    bound `p, q >= d` is now PROVED for all words, `hat S` included
    (item (35), FM-P35 (C)).
- **Stability (exact, checked).**  Take `d >= 1` and
  `kappa' = (|kappa| - kappa_1, kappa_2, ...)`, a boundary word with
  `kappa'_1 = |kappa'|/2`.  Then `m_(a,b)(h_kappa) = m_(a-2d,b)(h_kappa')`
  for all 75 partitions with `|kappa| <= 14` and `d >= 1`.
- **Consequence via the product formula (exact).**  With `s = r - d`,

      phi_r(h_kappa)/phi_r(empty) = Z_d(r) * sum_mu (-1)^b m_mu(h_kappa') [phi_s(V_mu)/phi_s(empty)] theta_mu(s),
      Z_d(r) = r(r-d) prod_(i<d) (r-i)^2,

  where

      theta_mu(s) = [(a+2d+2)/(a+2)] (s+2)^(q, rising) (s+3)^(p, rising) / [ (s+d+2)^(q+d, rising) (s+d+3)^(p+d, rising) ] > 0.

  - So every word with `d >= 1` is a positively reweighted copy of a
    boundary word at level `r - d`.
- **Observed (main agent, `reduce_d.py`).**  For all 45 `kappa` with
  `d >= 1` and `|kappa| <= 12`,
  `A_kappa(s) = R_kappa(s+d) / (Z_d(s+d) R_kappa'(s))` has no zero or
  pole in `(0, oo)` and is positive.  So the reweighting never flips the
  sign.
  - If that holds in general, FM3 and LP for the H-only sector reduce to
    words with `d <= 0` (interior and boundary).
  - `A_kappa` is not a product in general.  For example `(4,1,1)` gives
    `(s+2)(s+3)(3s^2 - 23s + 230) / (2(s+4)(s+5)^2(s+6)(s^2 - 3s + 12))`.

**(37) Transport formulation of the full cone (FM-STR, luna_max_jupiter;
main-agent exact checks `exactscan2.py`, `twolayer.py`).**
- **Setup.**  For a word `w`, put `F(mu) = (-1)^b W_r(mu) m_mu(w)` with the
  item (35) weights; FM3 at level `r` says `sum_mu F(mu) >= 0`.
- **H-only sector: grid transport.**  Negative cells can be matched to
  positive cells along the four grid neighbours `(a+-1, b+-1)`.
  - Exact integer max-flow: 86,376 (word, level) pairs are feasible.
    These are all H-only `|kappa| <= 24` (even) at `r <= 24`, plus hooks
    `(k,1^k)` for `k <= 30`.
  - The relative margin goes down to about `2e-15`, at the all-ones word
    `1^20`, `r = 24`.
  - No one-pass greedy rule works (item (35)).
- **`hat S` sectors: direct grid transport fails.**
  - Example: `h_1 hat S_3` at `r = 3` has cell values
    `F(0,0) = 1470`, `F(2,0) = 588`, `F(2,2) = -420`, `F(4,0) = 42`.
    The source `(2,2)` has no positive grid neighbour.
  - The move radius needed grows with the label; radius 7 is required at
    `kappa = (7)`, plus label 7, `r = 8`.
- **`hat S` sectors: two-layer transport.**
  - Expand
    `prod hat S_(p_i) = sum_gamma c_gamma h_gamma` using
    `hat S_p = 2h_p + 2h_(p-2) - h_(p-1) h_1`.
  - Layer `P` carries `A_mu = sum_(c > 0) c G_gamma(mu)`.  Layer `N`
    carries `-B_mu`, where `B_mu = sum_(c < 0) |c| G_gamma(mu)`.  Here
    `G_gamma` is the cell function of `h_(kappa cup gamma)`.
  - Allowed moves: the four grid neighbours inside each layer, and a
    vertical edge between `(P, mu)` and `(N, mu)`.
  - Exact max-flow is feasible in:
    - all 13,120 `hat S` pairs with total degree `<= 14`, `r <= 10`
      (FM-STR);
    - all 19,788 pairs over 15 plus patterns with total degree `<= 16`,
      `r <= 12` (main agent, independent code).
  - This is not implied by the H-only transport.  Layer `N` has net
    negative mass and must borrow across layers.
- **Conjecture T (transport form of FM3).**  For every integer `r >= 1`,
  every `kappa`, and every multiset of plus labels, the two-layer network
  (a single layer when there are no plus labels) admits a flow saturating
  all negative mass.
  - T implies FM3, and hence PCP and Q3.
  - T is a family of Hall inequalities: for every set `X` of negative
    nodes, the mass of `X` is at most the positive mass adjacent to `X`.
  - The Hall sets of smallest slack are staircase bands in the H-only
    data (FM-TR); FM-HALL is characterizing them.
- **Hall structure (FM-HALL, luna_max_saturn).**
  - **Proved (standard):** Hall's condition only needs source sets that
    are connected in the overlap graph, since slacks add over components.
  - **Band-minimizer conjecture.**  Order the odd cells by `p`
    increasing, `q` decreasing.  The minimum Hall slack is attained on a
    fixed interval ("staircase band") of this order; zero-demand cells are
    kept in the band.
    - Exact min-cut screen: all 7,860 entries with positive demand
      (H-only `|kappa| <= 16`, `r <= 16`, hooks `k <= 24`) agree.
    - Minimum slack 40, at `(2,1,1)`, `r = 2`, `S = {(1,0)}`.
    - The criterion is specific to FM3 arrays: it fails for arbitrary
      capacities already at `r = 4`.
  - The minimizing bands have sizes
    `1, 2, 4, 6, 9, 12, 16, 20, 25, ... = floor(P^2/4)`, i.e. full
    prefixes `{odd cells with p < P}`.
  - KILL: band kernels `R_(r,B)` are not positive combinations of ordinary
    level kernels.  Example: `r = 7`, `B = {(2,1), (3,2), (3,0)}` has a
    `V_(6,2)` term but no `V_(0,0)` term.
- **Single-cell Hall inequalities (FM-H1, luna_max_venus).**
  - Cell `(p,q) = (1,0)`, i.e. `mu = (1,1)`:

        m_11 <= (r+3)/(3(r-1)) m_00 + 2r/(3(r+2)) m_20 + 2(r-2)/(r+4) m_22.

  - Cell `(2,1)`, i.e. `mu = (3,1)`:

        m_31 <= 2(r+4)/(5(r-2)) m_20 + 6(r+2)/(5r) m_22 + 3(r-1)/(5(r+3)) m_40 + 9(r-3)/(5(r+5)) m_42.

  - Exact screens show no failure: the first inequality for even
    `|kappa| <= 24`, `r <= 60`, minimum slack 40; the second for
    `|kappa| <= 18`, `r <= 30`.
  - At `r = 2` the first inequality is exactly the open inequality T1
    (Corollary FM13), `3 I(k+(1,1)) <= 8 I(k) + 4 I(k+(2))`.  So even the
    first Hall cell contains the `r = 2` H-only core; no shortcut exists
    there.
- **Row-truncation positivity (main agent, `rowtrunc.py`).**  The prefix
  band `P` has slack

      sigma_P = sum_(p < P) F + sum_(p = P, p+q even) F,

  and `sigma_r` is FM3 itself.
  - `sigma_P >= 0` for every `P` in all 7,364 H-only pairs
    (`|kappa| <= 16`, `r <= 14`).
  - Also for all words with plus labels `(2)`, `(2,2)`, `(3,2)`.
  - It fails for some `hat S_3` words (8, first `(1,1,1)`, `r = 7`,
    `P = 3`) and `hat S_4` words (39).  Those sectors need the two-layer
    structure.
  - Plain truncation `sum_(p <= P) F >= 0` fails widely (2,831 H-only
    pairs, first `(1^4)`, `r = 5`).
  - As `r -> oo` with `P` fixed, `sigma_P` tends to the pure multiplicity
    inequality
    `sum_(p < P) (-1)^b (a+2)(b+1)/2 m_mu + (even row P) >= 0`.
    For example, `P = 2` gives
    `m_00 - 3 m_11 + 2 m_20 + 6 m_22 + 3 m_40 >= 0`, with coefficients
    `C_(p,q) = (p+q+2)(p-q+1)/2`.  An earlier version printed `3 m_22`;
    corrected by FM-H3.
- **Positive-mixture identity (FM-H3, luna_max_saturn; proved).**
  - Put `f_q(r) = prod_(i<q) (r-i)/(r+2+i)` and
    `g_p(r) = prod_(j<p) (r-1-j)/(r+3+j)`.
  - Let `beta_v` and `alpha_u` be the telescoping differences of `g` and
    `f`.  They are nonnegative and each family sums to 1.
  - For every fixed band `B`,

        sigma_B(h_kappa, r) / W_r(0,0) = sum_(v,u) beta_v(r) alpha_u(r) R_(B;v,u)(kappa),

    where `R_(B;v,u)` is the `r`-free clipped sum
    `sum_(p <= v, q <= u, (p,q) in B) (-1)^(p+q) C_(p,q) m_(p+q, p-q)`.
  - So every finite-level band slack is a probability mixture of `r`-free
    clipped multiplicity sums.
  - Screen: 656,360 exact row-prefix slacks (`|kappa| <= 20`, `r <= 40`),
    none negative.  The normalized slack is not monotone in `r` (example
    `(5,1)`, `P = 4`, `r = 9, 10, 11`).
  - The clipped sums themselves can be negative (`(1,1)`: `1 - 3 = -2`),
    and the `r = oo` limit inequalities do not imply the finite ones for
    arbitrary multiplicity vectors (`V_(1,1) + V_(4,0)`).  The special
    structure of `h_kappa` is essential.
- **KILL: inductive decomposition of row bands (FM-H2,
  luna_max_jupiter).**
  - `Q_r^[1] = (W_r(0,0)/6) Q_1 + (W_r(2,0)/6) Q_1 h_2` for every `r`.
  - The first nontrivial band `Q_3^[2]` is not a nonnegative combination
    of `Q_s u` with `s in {1,2}` and `u` among the 140 characters `1`,
    `h_k`, `hat S_p`, and products of two with index sum `<= 12`.
  - Exact integer separating functional: value `-6300` on the band, and
    `>= 0` on all 280 columns.
  - So row-band positivity does not follow from FM3 at lower levels on
    low-degree larger words.
- **T1 via the `p = gl(4)/sp(4)` action (FM-INJ, luna_max_venus): no
  progress.**
  - The contraction `A: Hom(U,M) -> M^Sp4` and the commutator map
    `B: Hom(U,M) -> Hom(Sym^2 W, M)` are jointly injective in the tested
    cases.  That only gives `m_U <= m_1 + m_ad`, not T1.
  - Any triple built from `A`, `A` composed with the Euler operator (which
    acts by the scalar `|kappa|`), and `B` has rank at most 13 < 15 at
    `(1^4)`.
  - The degree-raising fusion injection behind `t G_11 <= G_20` fails in
    the standard cyclic fusion grading at `(1,1)`.
  - T1 slack is 0 at one-part words, trivially, since all the relevant
    multiplicities vanish.
- **T1 via slot-dependent maps (FM-INJ2, luna_max_venus): a proof
  candidate.**
  - For each tensor slot `i`, put
    `A_i(f) = sum_a X_a^(vee,i) f(u_a)`, a map `Hom(U,M) -> M^Sp4`, and
    `B_i(f)(u ^ v) = X_u^(i) f(v) - X_v^(i) f(u)`, a map
    `Hom(U,M) -> Hom(Sym^2 W, M)`.
  - The block map
    `T: Hom(U,M)^(+3) -> (M^Sp4)^(+5) + Hom(Sym^2 W, M)` has entries
    `T_(r,c) = sum_i i^(r+1+5c) A_i` for `r <= 4` and
    `T_(5,c) = sum_i i^(5+c) B_i`.
  - Its rank, certified modulo `1,000,003`, equals `3 m_U` in all seven
    tested cases, so it is injective there:

    | `kappa` | `(m_1, m_U, m_ad)` | rank |
    |---|---|---|
    | `(1^4)` | `(3,5,6)` | 15 |
    | `(2,1,1)` | `(1,2,3)` | 6 |
    | `(2,2)` | `(1,1,1)` | 3 |
    | `(2,1^4)` | `(6,14,20)` | 42 |
    | `(1^6)` | `(14,30,40)` | 90 |
    | `(3,1,1,1)` | `(1,3,6)` | 9 |
    | `(2,2,1,1)` | `(3,7,10)` | 21 |

  - `(3,1,1,1)` is near-tight (slack 2).
  - **If `T_kappa` is injective for every `kappa`, T1 follows**, and with
    it the H-only `r = 2` sector.
  - **Main-agent independent reproduction (`inj/slotT2.py` to
    `slotT5.py`).**
    - The code is new and works on weight-0 subspaces.  It asserts that
      each `A_i(f)` is `Sp(4)`-invariant.
    - The multiplicities and injectivity agree on all of venus's cases.
    - **The specific power coefficients FAIL at the near-tight word
      `(2,2,2)`**: rank 8 < `3 m_U` = 9.
    - With random coefficients, `(2,2,2)` is injective for all three
      seeds tried.  So the power choice was degenerate.
    - At `(2,2,2)` the joint rank of all `A_i` is only `m_U = 3`, so the
      `B_i` channel is essential.
  - **Surviving candidate (generic-coefficient slot injection).**  For
    generic coefficient tensors `a_(r,c,i)`, `b_(c,i)`, the map

        T(f_0, f_1, f_2) = ( sum_(c,i) a_(r,c,i) A_i f_c )_(r = 0..4)  (+)  sum_(c,i) b_(c,i) B_i f_c

    is injective.
    - Exact mod-`p` screen with random coefficients: injective for every
      even partition with `|kappa| <= 6`, and so far for all tested
      size-8 words, including `(3,3,1,1)` (rank 21).  The run is
      continuing.
    - A uniform proof would need a transversality or rank argument for
      the pair of joint slot maps `(A_1..A_s)` and `(B_1..B_s)`.
    - FM-INJ3 (luna_max_venus, two primes) finds the fixed-power `T` of
      full rank on 11 further words:
      - `(2,2,2,2)` (45/45), `(3,3,1,1)`, `(4,2,1,1)`, `(4,1^4)`,
        `(3,2,2,1)`, `(4,4)`, `(4,3,1)`;
      - the near-tight `(3,1,1,1)` and `(3,1^5)` (90/90).
      It confirms rank 8 < 9 at `(2,2,2)` modulo both primes.
    - The joint `A` map alone cannot be the injection, since
      `4 m_1 = 24 < 45` at `(2,2,2,2)`.
    - Main-agent random-coefficient scan: `(3,2,2,1)` (27/27),
      `(3,2,1,1,1)` (51/51), `(4,1^4)` (12/12) and `(4,2,1,1)` (9/9) are
      injective.  No failure so far.
- **Level-step transport (main agent, `levelstep.py`).**
  - The exact level step is

        phi_r(w) = 4 phi_(r-1)(w h_2) + 8 phi_(r-1)(w) - 3 phi_(r-1)(w h_1^2).

    It comes from `(x-y)^2 = 3 hat S_2 + 2 - 2h_2` and matches Corollary
    FM13.
  - The three-layer network at level `r-1` has layers `4 w h_2`, `8 w`,
    `-3 w h_1^2`, with grid moves in each layer and vertical moves between
    layers.  Exact max-flow finds it feasible in all 1,440 cases (H-only
    `|kappa| <= 12`, `r = 2..10`).
  - So FM3 at level `r` has a transport certificate one level down.
    Iterating bottoms out at `r = 1` in the `Sp(4)` invariant-count model.
    This reorganizes FM3; it is not a proof.

**(38) Generic slot injection: a linear lift of Conjecture T (main agent,
2026-09-28; scripts `inj/slotHWg.py`, `hwscang.py`, `slotR.py`, `slotS.py`,
`hallstats.py`).**
- **Statement (Conjecture GSI).**  Fix `r` and a word.  For every cell `mu`
  of `Q_r` put `H_mu = Hom_Sp4(V_mu, M)` and give it `|c_mu|` copies, where
  `c_mu` is the `Q_r` coefficient (cleared of `2r(2r+1)`).  Negative cells
  are sources, positive cells targets.
  - Edge maps are the slot maps
    `S_i^(mu->nu)(f) = [ V_nu -> U (x) V_mu -> M, u (x) v -> X_u^(i) f(v) ]`
    for each tensor slot `i`, defined when `V_nu` is in `U (x) V_mu` (the
    grid neighbours of item (37)).  Here `X_u^(i)` is the `p = gl4/sp4`
    action on slot `i`.
  - Every (source copy, target copy) block of `T` is an independent generic
    combination of the available edge maps.
  - GSI: `T` is injective.  Injectivity gives
    `sum_neg |c_mu| m_mu <= sum_pos c_nu m_nu`, which is FM3 for the word at
    level `r`.  At `r = 2` this is the T1 slot map of item (37).
  - Conjecture T is the dimension shadow of GSI (subspaces `0` or all of
    `H_mu`).  GSI is stronger: it needs the linear Hall condition
    `sum_mu |c_mu| dim X_mu <= sum_nu c_nu dim(sum_(i, mu) S_i^(mu->nu) X_mu)`
    for all subspaces `X_mu` of `H_mu`.
- **`hat S` words: two-layer lift.**  Expand each `hat S_p` as
  `2h_p + 2h_(p-2) - h_(p-1) h_1` as in item (37), one module per choice.
  - Slot maps act inside each module.
  - Vertical edges join modules that differ in one factor.  They use the
    `G`-maps multiplication and `omega`-contraction
    `Sym^(p-1) W (x) W -> Sym^p W, Sym^(p-2) W`, and in the other direction
    polarization and `omega`-insertion.  Equivariance is checked for
    `p <= 5`.
  - **Kill (naive vertical lift).**  One map per vertical edge is not
    enough.  At `r = 1`, `kappa = (1,1)`, label `2`, the rank is 2 < 3,
    because the coefficient `2` in `2h_p` gives two copies but only one
    natural map.
  - **Repair (U-exchange).**  Add the composites
    `sum_a X_(Y_a^vee)^(i) o V o Y_a^(t)`, which pass the `U`-quantum
    between a tail slot `t` of the changed factor and any slot `i`.  These
    are `G`-maps from weight `mu` to weight `mu`.
- **Exact screens (mod 1,000,003; injectivity certified by full rank).**
  - `r = 2`, H-only (T1): EVERY even `kappa` with `kappa_1 <= |kappa|/2`,
    `|kappa| <= 12`, is injective, except `(1^12)` (out of memory range).
    - The size-12 scan completed on 2026-09-29, with 113 of 114 words
      done.  The last, `(2,1^10)`, reached rank 22770/22770 in 67,685 s
      at peak 47 GB.  `(1^12)` was skipped: its dense stage needs
      87.5 GB, over the 48 GB cap.
    - Size 14 (stopped 2026-09-30): 90 of the 105 words are injective,
      with no failure.  Nine were skipped at the 16.5 GB dense-stage cap:
      `(3,3,2,1^6)`, `(2^6,1,1)`, `(3,2^3,1^5)`, `(4,2,1^8)`,
      `(3,3,1^8)`, `(2^5,1^4)`, `(3,2,2,1^7)`, `(4,1^10)`, `(2^4,1^6)`
      (23.7--202 GB).  Six were not attempted, since each is larger:
      `(3,2,1^9)`, `(3,1^11)`, `(2,2,2,1^8)`, `(2,2,1^10)`, `(2,1^12)`,
      `(1^14)`.
  - `r = 3`, H-only: all 40 partitions of size `<= 8`, e.g. `(1^8)` at rank
    9100/9100.  The relative slack goes down to 9%, at `(2,1,1)` (80 vs 88).
  - `r = 4`, H-only: all 35 words run (size `<= 8`, the ones with many
    1s omitted), all injective.
    - Sources are `(1,1) x 378`, `(3,1) x 105`, `(3,3) x 35` and
      `(5,1) x 7`.
    - Examples: `(2^4)` at 9604/9604; `(3,2,1,1,1)` at 11319/11319;
      `(2,2,2,1,1)` at 17675/17675.  The last has slack 59, i.e. 0.3%.
  - **Structural constraint (main agent).**  At `r = 2` a slot-built
    `G`-map `Hom(U,M)^3 -> Hom(1,M)^5 + Hom(ad,M)` on the full `Hom`
    spaces goes between spaces of equal dimension, since
    `dim Q_2 = 5 - 15 + 10 = 0`.
    - Injectivity on the full space would make it a `G`-isomorphism and
      force `3 m_U = 5 m_1 + m_ad`.
    - T1 also fails for `V_k (x) N` with `N` arbitrary, e.g.
      `W (x) V_(2,1)` gives `-2`.
    - So no Lemma-JB-type (whole-module, one-slot-local) argument can prove
      GSI or T1.  A proof must use the invariant part, the fact that every
      factor is a symmetric power, and the rank-two (Pfaffian) structure.
  - `hat S` words, one label (U-exchange lift): labels `(2)`, `(3)`, `(4)`
    with `|kappa| <= 6` at `r = 1` are injective in all 90 cases (odd total
    degree is trivially zero).  `r = 2` passes on the cases run
    (`kappa = (), (1,1), (2), (2,2), (1^4)`, label 2; e.g. rank 224/224).
  - **KILL (two-label h-expansion lift).**  Labels `(2,2)` at `r = 1` fail.
    Every failing word below is short of `need` by the stated amount:

    | `kappa` | need | rank | short by |
    |---|---:|---:|---:|
    | `(2,1,1)` | 104 | 103 | 1 |
    | `(1^4)` | 216 | 213 | 3 |
    | `(2,2,2)` | 212 | 211 | 1 |
    | `(2,2,1,1)` | 404 | 400 | 4 |
    | `(2,1^4)` | 780 | 771 | 9 |
    | `(1^6)` | 1536 | 1515 | 21 |

    - Adding the `sp4` ("ad-exchange") composites
      `sum_b (Z_b^vee)^(i) o V o Z_b^(t)` does not change any of these ranks.
    - Diagnosis at `(2,1,1)`:
      - every edge map has full image rank;
      - multiplicity Hall holds with slack 4 (104 vs 108);
      - every proper subset of the four negative modules
        `(A,C), (C,A), (B,C), (C,B)` injects.
      So the obstruction is a single global linear relation, not a missing
      edge.  The likely source is that the two factors' vertical maps
      commute (double-complex structure of `hat S_2 (x) hat S_2`).
    - **The irreducible expansion fails the same way.**
      - The lift uses `hat S_p = V_(p,0) - V_(p-1,1) + V_(p-2,0)`:
        - `V_(p-1,1)` is realized as the joint kernel of multiplication and
          contraction on `Sym^(p-1) W (x) W`;
        - targets are projected by
          `1 - polar o mult / p - insert o contr / lambda_p`;
        - the vertical maps are the U-exchange `G`-maps.
      - One label passes (e.g. `(1^4)`, label 2: 5/5).
      - Labels `(2,2)` at `r = 1` fail:
        - `(2,1,1)`: 17 of 18;
        - `(1^4)`: 35 of 38;
        - `(2,2,2)`: 35 of 36.
      - These deficits (1, 3, 1) are the same as for the h-expansion lift.
        So the obstruction is intrinsic to vertical edges that change one
        factor at a time.
    - The `hat S` sector therefore needs maps between layers that differ
      in two or more factors.  The aggregated network of item (37) has
      such edges, and a `G`-map lift of them would need `ad`-type
      (two-quantum) exchange with `M`.
    - **RETRACTION (main agent, 2026-09-29; FM-CHK21 upheld).**  The Hodge
      and pair-map "repairs" below, and the evidence for Conjecture GSI-S,
      are int64-overflow artifacts.  They are withdrawn.
      - saturn's original FM-SLIFT script computed `proj @ mat @ F` with no
        mod-`p` reduction in between.  The intermediate entries reach about
        `3e19`, beyond the int64 range; the wrapped values act like random
        entries and inflate the rank.
      - Rerunning that exact script with the reduction inserted gives
        103/104 at `(2,1,1)`, labels `(2,2)`, for three seeds (104/104
        before the fix).
      - The saved FM-SLIFT5 map code (`cert/gsis_maps.py`, correctly
        reduced) run in the main harness also gives 103/104.
      - The local Hodge/pair maps are equivariant and exactly proportional
        between the two implementations (`cert/cmp_local.py`).  Only the
        assembly was wrong.
      - mercury's independent rebuild (FM-CHK21) found no repair in any of
        the six cases (103/104, 213/216, 49/50, 98/100, 62/64, 95/100),
        on two primes and two draws.
      - Still valid: the baseline-lift ranks computed with `slotS.py` (e.g.
        labels `(2,2)` at `r = 2` for `|kappa| <= 4`; labels `(2,2,2)` at
        `r = 1` up to `(2,1,1)`), and the deficits at `r = 1` with two
        labels.
      - So the two-label hat S lift is unrepaired.
      - **Recheck with the corrected harness, baseline lift only**
        (`cert/recheck.sh`; slot maps plus one-factor U-exchange vertical
        maps; no pair maps).
        - Injective:
          - `r = 2`: labels `(3,2)` for `kappa = (1), (3), (2,1), (1^3)`
            (e.g. 1665/1665), and labels `(3,3)` for `kappa = (), (2),
            (1,1)`;
          - `r = 1`: labels `(4,2)` and `(4,4)` for `kappa = (), (2),
            (1,1)`;
          - pure words (`kappa = ()`): labels `(2,2)`, `(3,3)`, `(4,2)`,
            `(4,4)`, `(5,3)` at `r = 1, 2, 3`, e.g. `r = 3`, `(3,3)`:
            1392/1392.
        - Not injective, all at `r = 1`:
          - labels `(3,2)`: `(2,1)` 49/50, `(1^3)` 98/100;
          - labels `(3,3)`: `(3,1)` 62/64, `(2,2)` 95/100;
          - labels `(2,2)`: `(1^4)` 213/216.
        - So the genuine hat S obstruction found so far sits at `r = 1`
          with two labels and nonempty `kappa`.
        - Rerunning with saturn's pair and `Gamma` maps added (correctly
          reduced) gives identical ranks in every case, including the five
          `r = 1` deficits.  These maps add no rank in any case tested;
          every pass is already a baseline pass.
      - **Full baseline screen (`inj/sscan.sh`, correctly reduced; all `kappa`
        with `|kappa| <= 6`).**
        - One label, `(2)`, `(3)` or `(4)`, at `r = 1, 2, 3`: injective in
          every case (about 30 words per label and level; odd total degree
          is trivially zero).
        - Two labels at `r = 2`, `(2,2)` and `(3,2)`: injective in every
          case.
        - Two labels at `r = 1`: the only failures.
          - Labels `(2,2)`: `(2,1,1)` 103/104, `(1^4)` 213/216, `(2,2,2)`
            211/212, `(2,2,1,1)` 400/404, `(2,1^4)` 771/780, `(1^6)`
            1515/1536.
          - Labels `(3,2)`: `(2,1)` 49/50, `(1^3)` 98/100, `(3,2)` 67/68,
            `(3,1,1)` 123/126, `(2,2,1)` 193/200, `(2,1^3)` 368/382,
            `(1^5)` 711/738.
        - Two labels at `r = 3`: not screened.  The run was stopped at
          75 GB resident (it was forcing the size-12 and size-14 T1 scans
          into swap) before any case was reported.  Single label `(4)` at
          `r = 3`: every `|kappa| <= 6` except `(1^6)` is injective;
          `(1^6)` was not run.
        - So the hat S obstruction to the slot-map lift sits exactly at
          level 1 with two or more labels.  At `r = 1`, FM3 is an
          invariant count, so this sector may need a separate argument.
      - **Level-1 two-label sector: exact reduction and a partial
        certificate (main agent, `xmn/xpq_check.py`, `conetest*.py`;
        check FM-CHK22 pending).**
        - Identity (sympy, all `p <= 7`):
          `hat S_p hat S_q = sum_(m = p-q, p-q+2, ..., p+q) hat S_m + X_pq`,
          with `X_pq = U_p(x)U_q(y) + U_q(x)U_p(y)
          = V_(p,q) + V_(p-2,q) - V_(p-1,q+1) - V_(p-1,q-1)`.  Here
          `chi_(a,b) = (U_(a+1)(x)U_b(y) - U_b(x)U_(a+1)(y))/(x-y)` also
          for non-dominant labels.
        - So `phi_1(hat S_p hat S_q h_kappa) = sum_m B_m + D_pq`, where
          `B_m = m_(m,0) + m_(m-2,0) - m_(m-1,1) >= 0` (FM53) and `D_pq`
          is the diamond `m_(p,q) + m_(p-2,q) - m_(p-1,q+1) - m_(p-1,q-1)`
          of `h_kappa`.  For a `GL(4)`-irreducible the `Sp(4)`
          constituents fill a rectangle in `(P,Q)` coordinates (FM51), so
          the diamond is a mixed second difference of a rectangle
          indicator.  Corrected by FM-SEC2 (saturn): it takes values in
          `{-1, 0, 1}` for `p - q >= 2`, is a first difference for
          `p - q = 1`, and equals `2(R(q,q) - R(q-1,q-1))`, with values
          in `{-2, 0, 2}`, on the diagonal `p = q`.
        - FM-SEC2 exact screen: 7,588 cases (`|kappa| <= 12`,
          `2 <= q <= p <= 8`), no negative total (3,761 zeros).  The
          per-shape payment fails first at `(p,q) = (2,2)`,
          `lambda = (2,1,1)` (combined coefficient `-1`).  A second
          per-shape test (expand one label as `2h_q + 2h_(q-2) -
          h_(q-1)h_1` and require a nonnegative coefficient on the FM51
          support of the other) fails badly, e.g. `-20` on `(4,2,2,2)`
          at `kappa = (1^8)` (main, `xmn/kostka_c.py`).
        - **KILL (per-`(m,n)` positivity).**  The `SU(2) x SU(2)`
          multiplicity array of `(x-y)^(2r) h_kappa` has a negative entry
          for every `kappa` with `|kappa| <= 8` at `r = 1, 2, 3` (67 of
          67 each; e.g. `(1,1)`: `-2` at `kappa = ()`).  So the diamonds
          cannot be made positive one at a time; they need the `hat S_m`
          terms.
        - **FM53-cone certificate.**  By LP, with supports forced to add
          because every generator restricts to a genuine
          `SU(2) x SU(2)`-module:
          - `hat S_p hat S_q` itself is not a nonnegative combination of
            `hat S_m (x) S_lambda W` and `S_lambda W` for any
            `p, q >= 2` (all `p + q <= 8`).
          - `hat S_2 hat S_2 h_kappa` is such a combination (hence
            `>= 0` by FM53) for every `kappa` with a part 1, and for
            `(3,2)`, `(3,3)`, `(5,2)`, `(4,3)`, `(5,3)`, `(4,4)`,
            `(2,2,2)`, `(4,2,2)`, `(3,3,2)`, `(2^4)`, `(3,2,2)`.
          - It is not for `(k)`, `k = 2..8`, nor for `(2,2)`, `(4,2)`,
            `(6,2)` (all partitions with parts `>= 2`, size `<= 8`).
          - The certified set is closed under adding parts, but its
            minimal elements look infinite (two-part words), so this
            gives no uniform proof by itself.
        - **Lemma BP (big part; PROVED by the main agent, formula checked
          in `xmn/bigpart_proof_check.py`; ACCEPTED by FM-CHK-BP,
          luna_max_mars).**  Let `p >= q >= 0` and
          `k >= p + q - 1`.  Then `X_pq (x) Sym^k W` is genuine, with
          `[X_pq (x) V_k : V_(a,b)] = mu(a-k, b) - mu(a-k, b+2)`.  Here
          `mu(i,j) = [|i| <= p, |j| <= q, i = p, j = q mod 2] +
          [|i| <= q, |j| <= p, i = q, j = p mod 2]` is the torus weight
          multiplicity of `X_pq`.  For `q = 0`, `X_p0 = hat S_p`.
          - Proof.
            - The character `U_p(x)U_q(y) + U_q(x)U_p(y)` of `X_pq` has
              torus weights `omega = (i,j)` with multiplicity
              `mu(i,j) >= 0`.  Brauer--Klimyk for the irreducible factor
              `V_(k,0)` gives `X_pq (x) V_k = sum_omega mu(omega)
              chi_((k,0)+omega)`, with the alternant convention.
            - Put `(A,B) = (k,0) + omega + rho`, `rho = (2,1)`.  For the
              first rectangle, `A >= k - p + 2 >= q + 1 >= |B|`.  For the
              second, `A >= k - q + 2 >= p + 1 >= |B|`.  So every point
              lies in `A >= |B|`, `A > 0`.
            - Points on the walls `B = 0`, `A = B` or `A = -B` contribute
              zero.  Points with `B < 0` reflect by `B -> -B` with sign
              `-1`, and no other Weyl element reaches the open chamber.
              Hence the multiplicity of `V_(a,b)` is
              `mu(a-k, b) - mu(a-k, -b-2) = mu(a-k, b) - mu(a-k, b+2)`.
            - For fixed `i`, `mu(i, .)` is a sum of indicators of
              symmetric intervals on one parity class, so it does not
              increase with `|j|`.  The difference is `>= 0`.
          - Check: the closed formula matches the direct decomposition
            on 115,692 entries (`p <= 8`, `0 <= q <= p`,
            `p+q-1 <= k <= p+q+7`), with no mismatch.
          - **Corollary (all labels, `r = 1`).**  Take labels
            `p_1, ..., p_l` and any `kappa` with some part
            `k >= p_1 + ... + p_l - 1`.  Then
            `phi_1(hat S_(p_1) ... hat S_(p_l) h_kappa) >= 0`.
            - `prod_j (U_(p_j)(x) + U_(p_j)(y)) = sum_(m,n) c_mn
              U_m(x) U_n(y)`, with `c_mn = c_nm >= 0` and
              `m + n <= sum p_j`.  So the product is a nonnegative
              combination of the `X_mn`.
            - Each `X_mn (x) V_k` is genuine by Lemma BP, so the
              product times `h_kappa` is genuine and has
              `phi_1 = dim(invariants) >= 0`.
            - FM53 is not needed here.
          - So the level-1 part of FM3 (any number of `hat S` labels)
            reduces to the small-part region: every part
            `<= sum p_j - 2`.
          - Equivalent `SU(2) x SU(2)` form (checked on 54 cases,
            `xmn/su2sq.py`).
            - Write `n_ab = [F : V_a (x) V_b]` for the restriction `F` of
              the word.  Then `phi_1 = n_00 + n_20 - n_11`.
            - Reason: `(x-y)^2 = sum_k (-1)^k Lambda^k p`, with
              `p = C^2 (x) C^2` the isotropy module of `Sp(4)/Sp(2)^2`.
            - So level-1 FM3 is the Euler characteristic of
              `Hom_H(Lambda^. p, F)` being `>= 0`.
            - With `hat S` factors, `F` is not a `g`-module, so no
              Chevalley--Eilenberg differential is available.
        - **Big-part screen (`xmn/xgen.py`, `xgen2.py`; FM-BIG1,
          luna_max_mars).**  `X_pq (x) Sym^k W` is a genuine
          representation whenever `k >= p + q - 1`.
          - Main agent: all `1 <= q <= p <= 7`, `k <= 15`.  Mars: an
            exact screen for `p <= 12` with no negative multiplicity above
            the threshold.
          - Below the threshold it fails in every tested case except
            `(p,q,k) = (3,1,1)`, where `X_31 (x) W = V_41` (mars).  So the
            statement is one-directional.
          - Mars's explicit eight-term alternation of the stencil
            `F_(p,q,k)` is in `agents/mars48.final.md`.  A proof for all
            parameters is still open.
          - The `q = 0` analogue, `hat S_p (x) Sym^k W` genuine for
            `k >= p - 1`, holds on mars's screen `p <= 8`.  Below it, the
            coefficient of `V_(p-1,k+1)` is `-1`.
          - Consequence.  If some part of `kappa` is `>= p + q - 1`, then
            `X_pq h_kappa` is genuine, and FM53 gives `B_m >= 0`.  So FM3
            at `r = 1` holds for `hat S_p hat S_q h_kappa`.
          - Conversely, among `|kappa| <= 10`, `X_pq h_kappa` fails to be
            genuine only when every part is `<= p + q - 2`.  For `(2,2)`
            these are 17 words: `(1^n)` for `n <= 8`, and words with at
            most two 2s.  For `(3,2)` there are 48; for `(3,3)`, 90; for
            `(4,2)`, 77.
          - The same region is where the two-layer slot lift fails: every
            baseline failure above has all parts `<= p + q - 2`.
          - Proof route: Brauer--Klimyk with the pyramid weight function
            `max(0, (k - |i| - |j|)/2 + 1)` of `Sym^k W`.  The diamond
            stencil is a difference of two second differences, so it
            vanishes away from the axes and the boundary.  The lemma
            becomes a finite piecewise-linear case check; not yet written
            out.
        - **Case `(p,q) = (2,2)` (curried slot map, `xmn/Ek.py`).**
          Here `X_22 = 2(V_(2,2) - V_(1,1))`.
          - The curried map `E_k: Hom(U, V_k) -> Hom(V_(2,2), V_k)`,
            built as in Lemma JB from `V_(2,2) = Sym^2_0 U`, has rank
            16/20, 45/50, 100/100, 175/175, 280/280, 420/420, 600/600 for
            `k = 1..7`.
            - At `k = 1` it kills the `V_(1,0)` summand of `U (x) W`;
              at `k = 2` it kills the `V_(1,1)` summand.
            - For `k >= 3` it is injective (checked to `k = 7`; the
              uniform proof is the `(2,2)` case of the big-part lemma).
            - The JB map was recomputed alongside, and its ranks agree
              with FM-CHK20.
          - So `m_(1,1) <= m_(2,2)` for every product with a factor
            `Sym^k W`, `k >= 3`.
          - Together with the FM53-cone certificates for `h_1` and
            `h_2^3`, and direct values at `()`, `(2)`, `(2,2)`
            (`D = 0`), the level-1 sector with labels `(2,2)` reduces to
            two exact rational LP identities plus the big-part lemma.
          - The two identities, exact over `Q` (`xmn/exactcert.py`; both
            sides as `SU(2) x SU(2)`-modules, where restriction is
            injective):
            - `hat S_2 hat S_2 h_1 = (1/10) hat S_1 s_31 + (1/5) hat S_1
              s_211 + (9/10) hat S_2 s_3 + (3/5) hat S_3 + (1/10) hat S_5`;
            - `hat S_2 hat S_2 h_2^3` is a combination of 20 terms
              `hat S_m s_lambda` with coefficients in `{5/7, 34/21, ...}`,
              all positive.
      - **Net-cell GSI: one lift for every sector (main agent,
        `inj/slotP.py`, `slotP2.py`, `slotPL.py`; FM-CHK23,
        luna_max_mercury, independent code: all 13 listed ranks
        reproduced at two primes and two coefficient draws, including
        the three net decompositions).**
        - Replace `Q_r` in Conjecture GSI by the net `Sp(4)`
          decomposition of the test representation
          `hat S_(p_1) ... hat S_(p_k) (x) Q_r`.  Slot maps act only on
          the `h_kappa` factors, and the `hat S` factors live in the cell
          weights.  Injectivity gives FM3 for the word.  With no labels
          this is GSI itself.
        - Example cells: `hat S_2 hat S_2 = 3V_00 - 3V_11 + 2V_20 +
          2V_22 - V_31 + V_40`; `hat S_3 hat S_3 = 3V_00 - V_11 + 2V_20 -
          2V_22 - V_31 + 2V_33 + 2V_40 - V_51 + V_60`.
        - Screen (all `|kappa| <= 6`, mod 1,000,003):
          - distance-1 edges (plain slot maps): injective for labels
            `(2,2)`, `(3,2)`, `(4,2)`, `(5,2)`, `(6,2)`, `(2,2,2)`,
            `(3,2,2)` at `r = 1`, and `(4,4)` at `r = 2`.  This includes
            every word where the two-layer lift failed.
          - distance-1 fails for `(3,1)`, `(5,1)`, `(4,3)`, `(5,3)`,
            `(3,3)` at `r = 1` and for `(3,3)` at `r = 2`.  Cause: for
            odd `q` the diamond `X_pq` has the opposite checkerboard sign
            to the `hat S_m` cells.  So a negative cell such as `V_22` in
            `hat S_3 hat S_3` has too few positive neighbours (at
            `(2,2)`: `m_22 = 1`, `m_33 = 0`).
          - distance-2 edges (composites of two slot maps): injective in
            every one of those cases.
          - distance-2 fails for `(4,4)` at `r = 1`: `(1^6)` 70/80,
            `(2,1^4)` 38/44.  The negative cell `V_33` lies three grid
            steps from every positive cell except `V_44`.
          - In general the diamond of `X_pq` sits `q - 1` steps from the
            `hat S_m` cells.  So with edges of bounded length, the
            small-part region (all parts `<= p + q - 2`) at `r = 1` is
            not covered for large `q`.
      - Audit of the main agent's code for the same risk:
        - in `slotR.py`, `c*(op@F1)` peaked at about `2e18`, below
          `9.2e18`, and old vs fixed images agree entry for entry
          (`cert/ovf_check.py`);
        - it is now reduced first;
        - `slotHWg.py` and `slotS.py` reduce every product.
        So the H-only GSI screens are unaffected.
    - **Partial repair (FM-SLIFT, luna_max_saturn) [WITHDRAWN, see
      retraction above].**
      - Use `A = V_20`, identified with `Lambda^2 U` for the orthogonal
        5-space `U`, and the `Sp4`-map `phi: A (x) U -> A`,
        `phi(w, u) = *(u ^ w)` (Hodge star).  It gives natural maps that
        change both factors: `(A,C) -> (B,A)` and `(C,A) -> (A,B)`.
      - With these added, the `(2,2)` h-expansion lift at `r = 1` is
        injective for every even `kappa` with `|kappa| <= 6`, up to
        `(1^6)` at 1536/1536.
      - The irreducible-expansion lift with the same two maps is still
        deficient at `(1^6)`: 265/270.
      - The baseline kernel at `(2,1,1)` is a circuit through all four
        source blocks; its coefficient-free form is not identified.
      - Open: labels `(3,2)`, `(2,2,2)`, `r >= 2`, all `kappa`.  The analogue
        for label `p` would use `V_(p,0) (x) U -> V_(p,0)`, which exists
        for all `p >= 1`.
    - **FM-SLIFT2 (saturn): mixed.**
      - The baseline lift (slot maps plus one-factor U-exchange vertical
        maps) passes:
        - labels `(2,2)` at `r = 2` for every even `|kappa| <= 4`, up to
          `(1^4)` at 3534/3534;
        - labels `(2,2,2)` at `r = 1` for `|kappa| <= 4` except `(1^4)`,
          which was not finished; e.g. `(2,1,1)` at 1824/1824.
      - It FAILS at `r = 1` for labels `(3,2)` and `(3,3)`:
        - `(3,2)`: `(2,1)` 49/50, `(1^3)` 98/100;
        - `(3,3)`: `(3,1)` 62/64, `(2,2)` 95/100, `(2,1,1)` 178/188,
          `(1^4)` 342/360.
        As at `(2,2)`, every proper subset of source nodes injects, so each
        deficit is a global circuit.
      - A general shared-constituent two-tail channel family was built, with
        equivariance checked.  It did not repair these, and it did not
        reproduce the accepted Hodge-star rank 104/104 at `(2,1,1)`, labels
        `(2,2)` (it gave 103).  It is not yet calibrated.
      - The failures concentrate at `r = 1` with two or more labels.
    - **FM-SLIFT3 (saturn).**
      - Calibration fixed.  The discrepancy was a coordinate-ordering bug
        in embedding the Hodge channel (`symbasis(1)` order).  After the
        fix, the accepted 104/104 and 216/216 are reproduced.
      - At labels `(3,2)`, `kappa = (2,1)`, `r = 1`, three independent
        coefficient draws each give rank 49/50 with different kernel
        lines.  So the deficit is a generic rank drop of this edge set,
        not one fixed relation.
      - Next candidate class (not yet built):
        `Gamma_(i,j,Phi) = sum_a rho_M(Y^a)^(i) o (1 (x) Phi) o
        rho_R(Y_a)^(j)`, i.e. a two-tail map `Phi` coupled to U-exchange
        with a slot of `M`.
    - **FM-SLIFT4 (saturn) [WITHDRAWN: overflow artifact, see retraction]: claimed repair on every screened case.**
      - Pair maps:
        - `Phi_32: A_3 (x) C_2 -> B_3 (x) A_2`, the composite
          `A_3 (x) W (x) W -> A_3 (x) U` (by `P_U`)
          `-> A_3` (Sym^3 action) `-> W (x) Sym^2 W` (polarization);
        - `Phi_33`: flip composed with the contraction `C_3 -> B_3`.
        Both pass the simple-root equivariance checks.
      - Added the `Gamma_(i,j,Phi)` class.  All pass at `r = 1`:
        - labels `(3,2)`, `|kappa| <= 3`, e.g. `(2,1)` 50/50 and `(1^3)`
          100/100;
        - labels `(3,3)`, `|kappa| <= 4`, e.g. `(3,1)` 64/64, `(2,2)`
          100/100, `(1^4)` 360/360;
        - the `(2,2)` regression for `|kappa| <= 6`.
      - **Conjecture GSI-S** (saturn) [its evidence is WITHDRAWN; the statement stands only as a proposal].  Take node spaces
        `H_(eps,mu) = Hom(V_mu, M_kappa (x) (x)_j R_(eps_j)(p_j))`,
        where `R_A = Sym^p`, `R_B = Sym^(p-2)` and
        `R_C = Sym^(p-1) (x) W`.
        - Copies are `|2^(#A+#B) Q_r(mu)|`, signed by the node mass.
        - Edges:
          - slot maps;
          - one-tail vertical maps (multiplication, contraction,
            polarization, insertion) and their U-exchange composites;
          - all two-tail shared-constituent maps `Phi`;
          - the `Gamma` exchanges.
        - Claim: a generic assembly has full source rank for every
          `kappa`, labels and `r`.
      - Pending independent check (FM-CHK21).  The ranks above come from
        inline code that was not saved.
    - **FM-SLIFT5 (saturn): GSI-S stress, two labels [pair-map ranks SUSPECT: to recheck with the corrected harness].**
      - Pair maps:
        - general `Phi^(AC)_(p,q) = Delta o mu o (1 (x) c_q)` and
          `Phi^(CA)_(p,q)`, from symmetric multiplication and
          comultiplication with symplectic contraction;
        - contraction-then-flip for equal labels;
        - the Hodge map at `p = 2`;
        - the `Gamma` exchanges.
      - All 35 two-label cases requested have full source rank:
        - `r = 2`, labels `(2,2)`: `|kappa| <= 4`, up to `(1^4)` at
          3534/3534;
        - `r = 2`, labels `(3,2)`: `|kappa| <= 3`, up to `(1^3)` at
          1665/1665;
        - `r = 2`, labels `(3,3)`: `|kappa| <= 2`;
        - `r = 1`, labels `(4,2)`, `(4,3)`, `(4,4)`: `|kappa| <= 2`.
      - Labels `(2,2,2)` at `(1^4)` did not finish.
      - The construction code is saved as `scratchpad/cert/gsis_maps.py`.
- **Lemma JB (PROVED; FM-INJ4, luna_max_venus; main-agent check).**  For
  every slot `i` with `kappa_i = k > 0`, the single map
  `B_i: Hom_Sp4(U, M) -> Hom_Sp4(Sym^2 W, M)` is injective.
  - Proof.
    - Currying gives `B_i = D_k (x) id_N`, with
      `D_k: Hom(U, V_k) -> Hom(Lambda^2 U, V_k)` and
      `D_k(t)(u ^ v) = X_u t(v) - X_v t(u)`.
    - `D_k` is `Sp4`-equivariant, and
      `U (x) V_(k,0) = V_(k+1,1) + V_(k,0) + V_(k-1,1)` for `k >= 2`, and
      `= V_(2,1) + V_(1,0)` for `k = 1`, is multiplicity free.  So `D_k` is
      injective iff it is nonzero on each summand.
    - Normalization (FM-CHK20 repair): take `U -> p` isometric,
      `X_(v^w)(z) = [omega(v,z) w - omega(w,z) v]/sqrt 2`, with
      Hilbert--Schmidt metrics.
      - The cross-Casimir `K = (C_(U (x) V_k) - C_U - C_(V_k))/2` has
        eigenvalues `k/2`, `-2`, `-(k+4)/2` on the three summands.
      - `D_k^* D_k = K(K+1)` then has eigenvalues `k(k+2)/4`, `2` and
        `(k+2)(k+4)/4`.
      - Without the factor `1/sqrt 2` they double; injectivity is
        unaffected.
    - On the `V_(k,0)` summand `{(X_a v)_a}`, `D_k` sends `v` to
      `([X_a, X_b] v)`.  This is nonzero for `v != 0`, because
      `[p,p] = sp4` and `V_k` is a nontrivial irreducible module.
    - `D_k (x) id_N` is injective on all of `Hom(U, V_k) (x) N`, hence on
      the `Sp4`-invariant part.
  - Main-agent check: `rank D_k = 5 dim V_k` (injective) for
    `k = 1, ..., 10`, mod 1,000,003.  Independent check (FM-CHK20,
    luna_max_mercury): ranks 20, 50, 100, 175, 280, 420 for
    `k = 1..6` mod 1009.  ACCEPTED after the two repairs above.
  - **Corollary.**  `m_(1,1)(M) <= m_(2,0)(M)` for every nonempty product
    `M = (x) Sym^(kappa_i) W`.
    - This strengthens the `m_(1,1) <= m_(2,0) + m_(0,0)` of FM51 for these
      products.
    - It holds in all 183 scanned words.
  - It is not enough for T1 by itself.  `T1 - JB = 5 m_1 - 2 m_U` is
    negative at `(2,2,2)` (`5 - 6`), so the `A` channel and at least two
    slots must interact.
- **Single-slot specialization (main agent, `sparse01.py`, `cyc.py`).**
  Generic coefficients are not needed in the slot direction.
  - Let each block of `T` use ONE slot with a generic scalar.  Then random
    slot assignments are almost always injective (e.g. 255/300 at
    `(2,2,2)`, 114/300 at `(2^4)`, 300/300 at `(1^6)`).
  - The fixed **cyclic rule** `T_cyc` is injective on all 15 words tested,
    up to `(1^8)` (rank 630/630), including `(2,2,2)`, `(2^4)`, `(3^4)`:

        T_cyc(f_0,f_1,f_2) = ( sum_c lambda_(rc) A_((r+c) mod s) f_c )_(r=0..4)  (+)  sum_c mu_c B_(c mod s) f_c

    with generic scalars `lambda`, `mu`.
  - With all scalars `= 1` the cyclic rule loses rank on symmetric words
    (`(2,2,2)`: 7/9; `(3^4)`: 81/90).  The scalars must break the slot
    symmetry; the slot pattern can be fixed.
  - **Reduction.**  By Lemma JB each `B_c` is injective.  Let
    `beta_3 = dim(B_0 H + B_1 H + B_2 H)`.  Then `T_cyc` is injective if
    (i) `3 m_U - beta_3 <= 5 m_1` and (ii) the five `A`-channels are
    jointly injective on the kernel of `sum_c mu_c B_c`.
    - At `(2,2,2)`: `beta_3 = 5`, and `9 - 5 = 4 <= 5`.
    - Measured (`beta3.py`, 18 words):
      - `beta_3 = m_ad` for most words and for every choice of three
        slots;
      - the exceptions are `(2,2,2)` (5/6), `(5,1^5)` (12/15), `(3,1^5)`
        and `(2,2,2,1,1)` (1 short for some triples).
      - (i) holds on all 18 words.
    - But since `beta_3 ~ m_ad`, condition (i) is essentially T1 itself.
      The split relocates the difficulty: it asks that the `B`-images of
      three slots span nearly all of `Hom(ad, M)`.  It does not reduce it.
- **KILL: crystal one-letter slot move (FM-INJ5, luna_max_venus).**
  - Setup.
    - Model `M` by the one-row KN crystals `B(kappa_1) (x) ... (x)
      B(kappa_s)`.
    - Move from a highest element of weight `mu`: change one letter of one
      factor by a weight of `U`, then raise to the highest element of the
      new component.  Keep the move if its weight `nu` is a positive
      `Q_r`-neighbour.
    - Give nodes their copy numbers and test Hall by exact max-flow.
  - Results.  It passes `r = 3` for `|kappa| <= 6` and `r = 2` for
    `|kappa| <= 6`.  It FAILS at `r = 2` on
    - `(2,2,2)`: flow 8/9;
    - `(3,3,2)`: flow 8/9;
    - `(2^4)`: flow 40/45.
  - The local-projection variant (keep only the `B(k)` component of
    `B(U) (x) B(k)`) fails already at `(1,1)`.
  - **Symmetry analysis (FM-SYM, luna_max_venus; corrects an earlier
    main-agent guess).**
    - Exact `S_kappa`-characters of `H`, `I`, `J`.  E.g. at `(2,2,2)`:
      - `H = sgn + Std`;
      - `I = sgn`;
      - `J = 2 x 1 + 2 Std`.
    - With trivial actions on copies and channels, an `S_kappa`-equivariant
      `T` has
      `rank <= sum_rho dim(rho) min(3 h_rho, 5 i_rho + j_rho)`.
      - This equals 7 at `(2,2,2)`, which matches the unit-scalar cyclic
        rank 7/9.
      - It is 9 at `(3,3,2)`, i.e. no deficit.  So slot symmetry does NOT
        explain the crystal failures (8/9 at `(2,2,2)` and `(3,3,2)`).
        Those are limitations of the one-letter move itself.
    - Scalars that depend on the slot alone still fail at `(2,2,2)` (8/9).
- **Explicit non-generic map (FM-SYM, venus; main-agent screen
  `expl.py`).**
  - Take the cyclic slot pattern `i_(rc) = (r+c) mod s`, `j_c = c mod s`
    with power scalars `lambda_(rc) = (i_(rc)+1)^(r+1+5c)` and
    `mu_c = (j_c+1)^(5+c)`.
  - This `T_exp` is injective on every word tested:
    - all 26 words of size `<= 8` with `m_U > 0`, up to `(1^8)` at
      630/630;
    - all 59 words of size 10 and 12 with `kappa_1 <= |kappa|/2` and
      `dim M_(1,1) <= 12000`, up to `(2,2,2,1^4)` at 729/729.
    That is 85 words with no failure.
  - It is a fixed, deterministic map.  Unlike GSI it can be attacked
    directly: its matrix entries are explicit powers.
  - **FM-EXP (venus).**
    - Two slots, PROVED: `T_exp` is injective for every word `(a,b)` at
      `r = 2`.
      - `dim Hom(U, V_a (x) V_b) = delta_ab`, and the generator is
        `F_a = omega(x,y)^(a-1) tr(Z C(x,y))`.
      - `A_0 F_a = A_1 F_a = ((a+4)/2) omega^a`, which is nonzero.
      - The first three `A`-rows then have coefficient matrix
        `((a+4)/2) [[1,64,1],[4,1,4096],[1,256,1]]`, with determinant
        `-785664 (a+4)^3/8`, which is nonzero.
      - This is a template only: T1 for two factors is already known
        (FM54).
    - The level-3 deterministic scheme uses slot
      `i(t,c) = (t+c) mod s` and exponent
      `e(t,c) = 1 + binom(t+c+1, 2) + c`, which is distinct per block.
      - It is injective on all 19 words of size `<= 6` at `r = 3`, up to
        `(1^6)` at 1225/1225.
      - The obstruction to the two-slot argument appears first at
        `(2,1,1)`, where `m_(1,1) = m_(3,1) = 2`.
    - **FM-EXP2 (venus): three slots, reduced but not closed.**
      - For `kappa = (a,b,c)`, `H` is spanned by
        `F_ij = omega_12^x12 omega_13^x13 omega_23^x23 C_ij`, with exponents
        fixed by the degrees.  So `dim H <= 3`, and there is no Fierz
        relation.
      - `T_exp` is injective on all 216 triples with `a, b, c <= 6`
        (45 of them with `dim H = 2`, 42 with `dim H = 3`).
      - Remaining step: when `dim H = 3`, the five `A`-coordinates have
        rank 5.  The `B`-coordinate must then be injective on the
        4-dimensional kernel.  This is a polynomial identity in the
        exponents, not yet proved.
      - Generic GSI at `r = 3` holds on all 22 partitions of size 8.
      - The deterministic level-3 map (main agent, `EXPLICIT=1 slotR.py`,
        reproducing venus's size-6 ranks) is injective on all 22 words of
        size 8, up to `(1^8)` at 9100/9100.  It is also injective on the
        tightest size-10 words: `(2^4,1,1)` 6275/6275 (slack 0.64%),
        `(3,2,2,1,1,1)` 4155/4155, `(2^5)` 3400/3400.
    - **FM-EXP3 (venus): symbolic three-slot minors; partial.**
      - Write `a = alpha + beta`, `b = alpha + gamma`, `c = beta + gamma`.
        Then `m_U = [alpha>=1] + [beta>=1] + [gamma>=1]` and
        `m_ad = 1, 3, 6` for `m_U = 1, 2, 3`.
      - Cases `m_U = 1, 2`: PROVED (conditional on the normalization
        check below).  The minors are nonzero multiples of
        `((x+4)/2)^3`, or `(...)^2 P(...)` with `P` positive on the
        feasible cone.
      - Case `m_U = 3`: one explicit `9 x 9` minor
        `D = -559872 alpha (alpha+beta+gamma+4) R` changes sign over the
        positive reals.
        - On the slice `(beta,gamma) = (4,1)`, a second minor survives at
          the sign change (the gcd is `559872(alpha+9)`).
        - There is no integer failure for `a, b, c <= 6`.
        - The general proof (some maximal minor nonzero at every positive
          integer triple) is open.
      - **Normalization reconciled (FM-EXP4).**
        - The screened map has endpoint coefficient `5/2`, so venus's
          table `A_i(F) = (alpha+beta+4)/2 Omega` at incident slots is the
          screened one.  FM-DEG1's `1/2` is corrected above.
        - The `B`-coefficients were checked against the code at
          `(2,1,1)`.
        - So the `m_U = 1, 2` proofs apply to `T_exp`.
    - **Theorem (three slots; FM-EXP5, luna_max_venus; certificate
      verified by the main agent).**  `T_exp` is injective for every
      three-factor word `kappa = (a,b,c)`.
      - Proof of `m_U = 3`.  Take the three `9 x 9` minors `D_1, D_2, D_3`
        (all five `A`-rows plus the `q`-rows `(11,12,22,23)`,
        `(12,13,22,23)`, `(12,13,23,33)`).  They satisfy
        `N D_1 - N_1 D_2 - N_2 D_3 = 559872 (alpha+beta+gamma+4)
        Q(alpha-1, beta-1, gamma-1)`, where:
        - `N = 1452003864120871163299950013`,
          `N_1 = 59424888149555732881404736`,
          `N_2 = 59440530383872825440470080`;
        - `Q` has 51 terms, total degree 5, and every coefficient
          positive.
        So some minor is nonzero on `[1, oo)^3`.
      - `m_U = 1, 2` are the boundary minors above.
      - Main-agent checks:
        - the verifier (`scratchpad/cert/s3_certificate.py`) reruns the
          determinant identity, the positivity of every coefficient, and
          the `6^3` grid;
        - `crosscheck.py` compares the basis-independent ranks (A-part,
          B-part, full, `m_U`, `m_ad`) of venus's symbolic matrix with the
          screened code on 9 three-factor words, and all match.
      - Scope.  T1 for three factors was already known (fifth pass,
        few-factor rational generating functions).  This is the first
        uniform injectivity theorem for the slot-map mechanism, and the
        first with multiplicity `> 1`.
      - Next instance: four slots, where Fierz relations appear, starting
        with `(1^4)`.
    - **Theorem T1-5 (main agent, 2026-09-29; verifiers
      `character_ring_iter/verify_t1_four_factors.py` and
      `verify_t1_five_factors.py`).**  T1, that is FM3 at `r = 2` in the
      H-only sector, holds for every word with at most five factors.
      This goes through the dimension inequality directly, not through
      GSI.  FM-EXP7 (venus) supplied the four-factor multiplicity
      formulas that suggested it.
      - ACCEPTED by FM-CHK-T15 (luna_max_neptune, independent code).
        - It cites Goodman--Wallach Thm 5.2.2, Cor 5.2.5 and
          Thm 12.2.15 for the fundamental theorems, and
          Buchsbaum--Eisenbud (Amer. J. Math. 99 (1977)) Thm 2.1 for
          the resolution.
        - Series re-checked on 1,001 four-part and 3,003 five-part
          compositions (size `<= 10`).
        - Both certificates re-expanded exactly.
      - Step 1.  `T1 = 8 m0(kappa) + 4 m0(kappa,2) - 3 m0(kappa,1,1)`,
        where `m0` counts invariants and `W (x) W = 1 + U + Sym^2 W`.
      - Step 2 (fundamental theorems for `Sp(4)`).  The invariants of `n`
        vectors are `C[omega_ij]/I_6`, with `I_6` the ideal of `6 x 6`
        Pfaffians.
        - `I_6 = 0` for `n <= 5`, and `I_6 = (Pf)` for `n = 6`.
        - For `n = 7`, Buchsbaum--Eisenbud gives the resolution
          `0 -> R(-(2^7)) -> (+)_l R(-(1^7) - e_l) ->
          (+)_k R(-(1^7) + e_k) -> R`.
        - Extract `t^2` from `n = s+1` and `t1 t2` from `n = s+2`, with
          `Pi_s = prod_(i<j<=s) (1 - u_i u_j)^(-1)`.  This gives
          `sum_kappa T1(kappa) u^kappa = N_s Pi_s`, with:
          - `N_4 = 5 + p_2 - 2 e_2 + 3 e_4`;
          - `N_5 = 5 + p_2 - 2 e_2 + 3 e_4 - e_1 e_5`.
        - Checked against direct Weyl-character values of T1 on 922
          four-part and 1,292 five-part compositions (size `<= 14` and
          `<= 12`), with no mismatch.
      - Step 3 (positive decomposition, found by LP and verified
        exactly).  `N_s = sum_i c_i u^(m_i) prod_(e in F_i) (1 - u_e)`
        with `c_i > 0` and `F_i` sets of edges of `K_s`.  So
        `N_s Pi_s = sum_i c_i u^(m_i) prod_(e not in F_i)
        (1 - u_e)^(-1)` has nonnegative coefficients.
        - `s = 4`: nine terms, all with coefficient 1:
          `N_4 = w12 + u1^2 w34 + w13 w14 + u3^2 w12 w24 + u4^2 w13 w23
          + u2^2 w14 w34 + w12 w24 w34 + w13 w23 w24 + w14 w23 w34`,
          where `wij = 1 - u_i u_j`.
        - `s = 5`: 56 terms, `u`-degree `<= 8`, rational coefficients
          with denominators 7, 19, 133 (listed in the verifier).
          Degree 6 is infeasible.
      - Scope.
        - T1 for three factors was known.  Four and five factors are new,
          for all degrees at once.
        - The method is the same at every level `r` and in the `hat S`
          sectors: a rational series in the factor variables whose
          numerator needs a positive decomposition.  But the invariant
          rings on `s + 2(r-1)` or more vertices carry the full Pfaffian
          resolutions.
        - From `s = 6` on, the Pfaffian ideal enters `m0(kappa)` itself.
          A uniform-in-`s` decomposition is open.
    - **Theorem T1-6 (main agent; verifier `verify_t1_six_factors.py`;
      ACCEPTED by FM-CHK-T16, luna_max_neptune).**  T1 holds for every
      word with at most six factors.
      - The checker supplied the citation for the degree bound:
        Bruns--Herzog, *On the computation of a-invariants*, Cor. 1.7
        (`a = -rn` for Pfaffian ideals of size `2r+2`, `2r < n`), with
        the Gorenstein property from Kleppe--Laksov.
      - It re-derived `deg N_6 <= 26` and rechecked `N_6` on all 74,613
        six-part compositions of size `<= 16`.
      - `N_6 = 5 + p_2 - 2e_2 + 3e_4 - e_1e_5 - 4e_6 + e_2e_6 - e_6^2`.
        It comes from exact T1 values of all 6-part compositions through
        `u`-degree 26, with nothing above degree 12.
        - Degree 26 is the a priori bound: the 8-vertex Pfaffian ring is
          Gorenstein with `a = -2n`, so `K_8` has `x`-degree `<= 24`,
          and the `t_1 t_2` extraction adds at most 2.
        - The a-invariant formula matches `n = 6` (principal Pfaffian)
          and `n = 7` (Buchsbaum--Eisenbud top term `x^(2^7)`).
      - Cycle decomposition (exact):
        `N_6 = (1/12) [ sum over the 60 Hamiltonian cycles C of K_6 of
        prod_(e in C) (1 - u_e) + sum over vertices v and 5-cycles C
        on the other vertices of u_v^2 prod_(e in C) (1 - u_e) ]`.
      - So `12 T1(kappa) = sum_C G_(K_6 - C)(kappa) + sum_(v,C)
        G_(K_6 - C)(kappa - 2e_v)`.  Here `G_X(kappa)` counts multigraphs
        with edges in `X` and degree sequence `kappa`, so T1 is a count
        of pairs (cycle, multigraph avoiding its edges).
      - Five factors in the same form: `N_5 = (1/12) [ sum over
        Hamiltonian paths P of prod_P (1 - u_e) + sum over (v, 4-vertex
        path P off v) of u_v^2 prod_P (1 - u_e) ]`.  This comes from the
        symmetric LP, as an alternative to the 56-term certificate.
      - Seven factors: `N_7 = N_6 + e_7 (4e_1 - 3e_3 + 4e_5 - 5e_7)`.
        It comes from exact T1 values through degree 22 (no terms in
        degrees 16--22; the a priori bound is 38); it terminates
        at degree 14, and its absolute coefficient mass by degree
        (5, 49, 105, 175, 175, 105, 49, 5) is symmetric.  No positive
        decomposition of this form exists with `u`-degree `<= 14`
        (symmetric LP over 263,370 columns, infeasible).  Degree 16
        (675,806 columns) is also infeasible.
      - Direct T1 scans (Weyl character formula, `xmn/t1scan7.py`): all
        7-part partitions of size `<= 30` and 8-part partitions of size
        `<= 32`, no negative value.  The only zeros are the forced ones
        (`d >= 2`), and every `d = 1` word has `T1 = 1`.
      - **Observation (Schur form).**  In the Schur basis every `N_s`,
        `4 <= s <= 7`, uses only shapes with at most two columns.
        - `N_4`, `N_5`, `N_6` are exactly the truncations of `N_7` to
          4, 5, 6 rows.  So there is one numerator
          `N = sum c(a,b) s_(2^a 1^b)`, and `N_s` is its restriction to
          `s` variables.
        - Known coefficients:
          - `c(0,0) = 5`, `c(1,0) = 1`, `c(0,2) = -3`, `c(0,4) = 3`,
            `c(1,4) = -1`;
          - `c(2,4) = 1`, `c(0,6) = -5`, `c(6,0) = -1`;
          - `c(1,6) = 5`, `c(3,4) = -3`, `c(5,2) = 3`, `c(7,0) = -5`;
          - all others with `a + b <= 7` are 0.
        - These satisfy `c(7-a-b, b) = -c(a,b)`.
        - If the two-column form holds in general, `N_s` has degree
          `<= 2s`, and `N_7` has no terms above degree 14.
        - **Correction (FM-GEN1, venus).**  The two-column form fails
          at `s = 8`.
          - Exact T1 data through degree 16 give coefficient `-5` at
            `m_(3,1^7)` in `N_8`, and no Schur function with at most two
            columns has that monomial.
          - The degree-16 part of `N_8` also contains `55 m_(1^8)`,
            `55 m_(2^8)` and `-20 m_(2^2,1^6)`.
          - So the two-column form holds only through 7 rows.  Whether
            `N_7` has terms above degree 22 is decided only by the
            a priori bound (degree 38).
      - Diagnostic (main agent, `xmn/posdec7K.py`).
        - The 7-vertex Pfaffian K-polynomial
          `K_7 = 1 - e_6 + e_1 e_7 - e_7^2` IS product-certifiable.
          With all monomials at degree 14 the symmetric LP is feasible,
          and it uses odd, edge-like monomials such as `u_3 u_5`.  This
          is consistent with a Stanley decomposition from a shelling of
          the 3-nesting-free complex.
        - So the class is not too weak for Pfaffian rings as such.
          `N_7` lies outside it at degree `<= 16` for a reason of its
          own.
        - A product-class certificate is the same as a lift of `N_7` to
          edge variables `x_e` whose series over `prod (1 - x_e)` is
          nonnegative edge by edge.  That is a nonnegative weight on
          multigraphs summing to `T1(kappa)` over each degree class.
          For `s = 7`, no lift of bounded degree `<= 16` exists.
      - **Seven factors by slicing (main agent, `xmn/slice7.py`,
        `posdec6p1.py`, `slice_exact.py`; slices `k = 1..4` ACCEPTED by
        FM-CHK24, luna_max_mercury).**
        - Mercury re-derived the slice formula and the symmetry
          argument.
        - It re-solved the four certificates exactly from their
          supports: ranks 5, 29, 64, 130, zero residual, all
          coefficients positive.
        - The `u_7^k` coefficient of `N_7 Pi_7` is `N^(k) Pi_6`, with
          `N^(k) = sum_(j <= 2) [u_7^j] N_7 * h_(k-j)(u_1..u_6)`.
        - Each slice is a six-variable problem.  So the product class in
          six variables, with a certificate that may depend on `k`, is
          enough.
        - Exact `S_6`-symmetric certificates:
          - `k = 0`: this is `N_6`, the two cycle terms;
          - `k = 1`: 5 terms at degree 15, all on the 5-cycle (vertex
            factor `u_v`, `u_v^3`, `u_a u_b u_v^3`) or the Hamiltonian
            6-cycle (`u_v`, `u_a u_b u_v`);
          - `k = 2`: 29 terms at degree 18, rational, minimum `1/720`,
            reconstructed exactly.
          - `k = 3`: 64 terms at degree 19, rational, minimum
            `29/108180`, reconstructed exactly with flint `fmpq_mat`
            (`xmn/slice_exact_flint.py`; rank 64, no free parameter,
            zero residual).
          - `k = 4`: 130 terms at degree 20, reconstructed exactly
            (rank 130, no free parameter).
        - So T1 holds for every seven-factor word whose smallest part is
          `<= 4`, by symmetry of T1 in the parts.
        - FM-GEN3 (venus): `N^(k) = N_6 h_k + n_1 h_(k-1) + n_2 h_(k-2)`,
          where
          - `n_1 = -2E1 + 3E3 - E1E4 - 5E5 + E2E5 + 5E1E6 - 3E3E6 +
            2E5E6`;
          - `n_2 = 1 - E4 + E1E5 - E5^2 + E6(4 - 3E2 + 4E4 - 5E6)`;
          - `E_j = e_j(u_1..u_6)`.
          Venus re-derived the `k = 0, 1` certificates independently.
          No family valid for all `k` was found, and `N_6 h_k` alone
          fails at `k = 1` (its `u_6` coefficient is `-2`).
        - The degree needed above the slice's own degree is 0, 2, 4, 4,
          4 for `k = 0..4`.  `k = 3` is infeasible at degree 17 (phase-1
          slack 0.0027).  The certificate size grows: 2, 5, 29, 64,
          130 terms.
        - `k = 5` is infeasible at degree 21 (1,241,232 columns, phase-1
          slack 0.00015, 2.9 h).  So the gap exceeds 4 at `k = 5`, and
          each higher slice costs more.
        - As a route to all seven-factor words, slicing is exhausted:
          it proves one slice at a time and has no uniform family
          (FM-GEN3).
        - KILL (main agent, `xmn/shift7.py`): a shifted-region
          certificate.
          - For the complement region (all parts `>= 5`), the shifted
            series `sum_lambda T1(5*1 + lambda) u^lambda` times
            `prod (1 - u_i u_j)` does not terminate.  Its absolute
            coefficient mass grows from `8 * 10^4` at degree 1 to
            `7 * 10^8` at degree 15.
          - So there is no low-degree numerator to decompose.
          - T1 is large there: the minimum is 1,563 over all 316 sorted
            `lambda` with `|lambda| <= 15`.  This region needs a
            lower-bound argument instead.
          - For `k = 2` at degree 16 the slack is 0.00085.
        - HiGHS returns status 15 on the plain formulation; the phase-1
          slack formulation is reliable (the `k = 1` control has slack
          0).
      - **KILL (FM-COMB1, luna_max_venus): natural extensions of the
        cycle count to `s >= 7`.**
        - The count sums the 6-cycle and (vertex, 5-cycle) terms over
          all 6-subsets of `[s]`, with `H` any multigraph on `K_s`
          avoiding the cycle.  Normalization `c = 12` is fixed by
          `s = 6`.
        - It fails at `(2,1^6)`: 5,400 against `12 T1 = 480`.  At
          `(1^8)` it gives 67,200 against 1,080.
        - Requiring `H` 3-nesting-free, or `H + Gamma` 3-nesting-free,
          fails already at `(1^6)` (224 and 168 against 240).
        - The `s = 5` path control passes.  T1 itself is not falsified.
      - **Conjecture EM2 (smallest-edge descent; main agent,
        `xmn/mono.py`, `mono2.py`, `em2.py`).**
        - Statement: for every partition `kappa` with at least two
          positive parts, `kappa != (1,1)`,
          `T1(kappa) >= T1(kappa - e_(s-1) - e_s)`, where `kappa_(s-1)`,
          `kappa_s` are the two smallest positive parts.
        - Evidence: no failure on all 3,737 partitions with `<= 12` parts
          and size `<= 24`.  The minimum margin is 0, at `(3,1) -> (2)`.
          - Also no failure on all 10,316 partitions with exactly 7
            parts and size `<= 40`, and all 5,830 with exactly 8 parts
            and size `<= 36` (`xmn/em2big.py`).
          - The only zero margins there are the forced ones, e.g.
            `(10,1^6)` and `(11,1^7)`.
        - Other rules fail:
          - removing the edge between the two largest parts fails 79
            times in 444 words, first at `(2,1,1)`;
          - removing the edge between the largest and smallest parts
            fails on `(k+1,k,1)` for every `k`;
          - monotonicity under adding an arbitrary edge fails in 202 of
            1,980 cases.
        - Consequences.
          - EM2 plus the base cases (single parts, and `T1(1,1) = 3`)
            gives T1 for every word by descent.
          - Each step lowers the two smallest parts, so every word
            descends to one with fewer factors.  Given Theorem T1-6, EM2
            is needed only for words with `>= 7` positive parts.
        - Independent check (FM-CHK25, luna_max_mercury): ACCEPT on all
          four points.
          - EM2 reproduced on 6,694 inequalities, and the descent
            implication confirmed.
          - The `D_(a,b)` identity confirmed.
          - All six verifiers of commit `c963e28` rerun.
        - FM-EM2A (luna_max_mars), seven parts, all `>= 5`.
          - Odd total degree gives 0 = 0, and `d >= 1` follows from
            Lemma D1.  So only the balanced cone `d <= 0` is open.
          - Exact screen: 564 sorted tuples, degree `<= 50`, no negative
            margin.  The smallest positive margin is 812, at
            `(20,5^6)`.
        - FM-LIT6 (luna_max_neptune): no theorem in the literature
          implies EM2.
          - Sources checked: Jonsson--Welker, de Mier, Ghorpade--
            Krattenthaler, Fayers, White, Huh et al., Branden--Huh,
            Lam--Pylyavskyy, Enright--Willenbring.
          - Reformulation.  `R = C[omega]/Pf_6` is a domain, so
            multiplication by `omega_ij` is injective.  Each increment
            `delta = H(kappa) - H(kappa - e_i - e_j)` is the Hilbert
            function of `R/(omega_ij)`, and likewise on the extended
            vertex sets.
          - So EM2 is T1 for the quotient ring `R/(omega_ij)`, with
            `i, j` the two smallest-degree vertices:
            `3 delta_11 <= 8 delta_0 + 4 delta_2`.
          - In representation terms this is T1 for
            `h_kappa'' (x) (Sym^a W (x) Sym^b W)/omega = h_kappa''
            (x) D_(a,b)`.
        - FM-EM2D (luna_max_mars): cone generating function;
          conditional.
          - The exact cone transform is
            `C(z) = Omega_(>= 0) [ z_0^(-5) N_7(v) / prod_(e != 67)
            (1 - v_e) ]`, with `v` the cone substitution.
          - The Omega elimination was not completed: edge factors such
            as `v_1 v_3 = z_1 z_3 / z_2` have mixed-sign exponents.
          - EM2 margins are negative off the cone (for pairs other than
            the two smallest), so any certificate must be cone-restricted.
            This class is new; the product class cannot work.
        - FM-EM2B (luna_max_saturn), all `s`: conditional.
          - Paying cell by cell fails: at `a = b = 1`,
            `kappa' = (2)`, `V_(1,1) h_2` gives `-2` and `V_(2,0) h_2`
            gives `+3`.
          - The Kostka-summed inequality is open.
        - KILL (main agent, `xmn/em2conv.py`): diagonal convexity.
          - The telescope `h_a h_b = sum_(i=0..b) D_(a-i,b-i)` makes EM2
            say every step is nonnegative.
          - Convexity `Delta(kappa',a,b) >= Delta(kappa',a-1,b-1)` would
            reduce EM2 to `b = 1`.  It fails in 36 of 501 cases, first at
            `(2,2,2)` (`0 < 1`).
          - The `b = 1` case itself holds on all 1,438 words tested
            (size `<= 22`, `<= 9` parts).
        - KILL (FM-EM2C, luna_max_venus): per-constituent and per-Pieri
          proofs of the `b = 1` case.
          - The ordering hypothesis does not make each `GL(4)`
            constituent pair nonnegatively: at `a = 1`, `kappa' = (1^4)`,
            `s_(1^4) D_(1,1)` gives `-2`.
          - Termwise Pieri also fails: `s_111 h_1 = s_211 + s_1111` gives
            `5` and `-2`.
          - The aggregate `b = 1` inequality is open.
        - Equivalent form.  With `a >= b` the two smallest parts,
          `h_a h_b - h_(a-1) h_(b-1) = D_(a,b) = sum_(j=0..b)
          V_(a+b-j, j)`, a genuine module by Pieri.  So EM2 says
          `<h_kappa' (x) D_(a,b), Q_2> >= 0` whenever every part of
          `kappa'` is `>= a`.
          - The small factors enter only through the genuine row
            `D_(a,b)`.
        - KILL (main agent, `xmn/em2r.py`, `xmn/descent_r.py`): EM2 holds
          only at `r = 2`.
          - At level `r`, the descent `phi_r(kappa) >=
            phi_r(kappa - e_(s-1) - e_s)` fails for 14, 87 and 142 of the
            517 partitions of size `<= 16` at `r = 3, 4, 5`.
          - At `r = 3` every failure has smallest part 1, e.g. `(3,1)`
            (`-4`), `(2,1,1)` (`-6`), `(2,2,1,1)` (`-9`).  From `r = 4`
            large parts fail too: `(5,5,4,4)` at `r = 4` (`-17`),
            `(6,5,5)` at `r = 5` (`-28`).
          - Six other descents fail at every `r = 2..5`:
            - the edge between the two largest parts;
            - the edge between the largest and smallest parts;
            - merging the two smallest parts;
            - lowering the largest or the smallest part by 2;
            - dropping an equal smallest pair.
          - So EM2, and any quotient-ring lift of it, can close at most the
            `r = 2` H-only sector.  The level-uniform candidate is still
            GSI.
        - **Conjecture QGSI (quotient slot injection; main agent,
          `inj/qgsi.py`, `inj/qscan.py`).**  A linear lift of EM2.
          - Realize `h_kappa'' (x) D_(a,b)` as `ker Delta`, where `Delta`
            is the `omega`-contraction of the last two slots (parts `a`,
            `b`).  Slot maps `A_i`, `B_i` on the inner slots `i <= s-2`
            commute with `Delta`, so they act on this quotient.  Checked:
            `Delta` kills every `A_i` image.
          - QGSI: the generic `3 -> 5 + 1` combination of the inner slot
            maps is injective on `Hom(U, ker Delta)^3`.  QGSI implies EM2.
          - On the cone (`a, b` the two smallest parts), no failure with
            three or more inner slots.  Complete scan: all 67 cone words
            with 5 to 10 parts and size `<= 12` (`<= 9` parts at size 12;
            parts 5/6/7/8/9/10: 24/19/11/8/4/1 words) are injective.  The
            heaviest, `(2,2,2 | 1^6)`, took 4.3 h.
          - With two inner slots it fails: `(2,2 | 1,1)` 16/18,
            `(2,2 | 2,2)` 22/24, `(3,2 | 2,1)` 18/21, `(3,3 | 1,1)` 16/18,
            `(4,2 | 1,1)` 8/9.  `(1,1 | 1,1)` passes.  EM2 is needed only
            for `>= 7` parts, i.e. `>= 5` inner slots.
          - Off the cone it can fail even when the dimensions allow
            injectivity: `(1^4 | 3,1)` has margin 2 and rank 47/48, and
            `(2,1,1 | 2,2)` has margin 1 and rank 41/42.  So QGSI is not a
            generic maximal-rank statement.  `(1^6 | 3,1)`,
            `(1^6 | 2,2)` and `(1,1,2,2 | 2,2)` pass.
          - The ranks agree for three random seeds.
          - An earlier run that did not reduce two sparse products mod `p`
            overflowed `int64` and reported every word injective.  It is
            discarded (`qscan_overflow_suspect.*`).
          - Scope: `r = 2`, H-only (EM2 fails at `r >= 3`).  A proof of
            QGSI faces the GSI problem with fewer slot maps, so no new
            mechanism is in view.
      - **Lemma D1 (FM-GEN1, venus; proved).**  `T1(kappa) = 1` for
        every word with `kappa_1 = |kappa|/2 + 1`, for any number of
        factors.
        - At multidegrees `(kappa,2)` and `(kappa,1,1)` every multigraph
          is a star at vertex 1.  The Pfaffian ideal has no component
          there.
        - So `m_00 = m_11 = 0`, `m_20 = 1` and `T1 = 1`.
    - **Theorem S2-2 (main agent; verifier
      `character_ring_iter/verify_shat_pair_two_factors.py`; S2-2 and
      S2-3 ACCEPTED by FM-CHK-S2, luna_max_neptune, with independent
      numerator derivation and Weyl checks on 1,820 + 6,188
      compositions).**  FM3 at
      `r = 1` holds for every word `hat S_p hat S_q h_a h_b`, for all
      `p, q, a, b >= 0` at once.
      - The labels are treated as variables:
        `sum_p hat S_p t^p = (2 + 2t^2) H(t) - t H(t) h_1(z)`, where `z`
        is an extra degree-1 vertex.  Only invariant rings on `<= 6`
        vertices enter.
      - This gives
        `G(t1,t2,u1,u2) prod_(K_4)(1 - x_i x_j) = N`, with
        `N = (2+2t1^2)(2+2t2^2) - t1(2+2t2^2) h1 - t2(2+2t1^2) h1 +
        t1 t2 (h1^2 + 1 - e4)`.  The closed form matches exact
        `phi_1` values through degree 14.
      - `N` is a sum of ten terms `x^m prod_(e in F)(1 - x_e)`, all
        with coefficient 1, so `G` has nonnegative coefficients.
      - This covers the small-part region for these words.
    - **Theorem S2-3 (main agent; verifier
      `character_ring_iter/verify_shat_pair_three_factors.py`).**  FM3 at
      `r = 1` holds for every word `hat S_p hat S_q h_a h_b h_c`, for all
      `p, q, a, b, c >= 0`.
      - The numerator uses invariant rings on `<= 7` vertices, so
        Buchsbaum--Eisenbud enters:
        `N = (2+2t1^2)(2+2t2^2) - t1(2+2t2^2)(e1 - e5) -
        t2(2+2t1^2)(e1 - e5) + t1 t2 (e1^2 + 1 - e4 - e1 e5)`,
        with `e_k` in `(t1, t2, u1, u2, u3)`.  It matches exact `phi_1`
        values through degree 12.
      - It has a positive decomposition with 64 terms, rational
        coefficients (denominators dividing 45), and `x`-degree
        `<= 12`.  Degrees 8 and 10 are infeasible.
    - Level-1 two-label sector, current state: proved for words with at
      most three `h`-factors (Theorems S2-2, S2-3), and for words with a
      part `>= p + q - 1` (Lemma BP).  Open: four or more small parts.
      - Four `h`-factors, closed form (8 vertices, via the `N_6`
        extraction):
        `N = (2+2t1^2)(2+2t2^2)(1-e6) - (t1(2+2t2^2) + t2(2+2t1^2))
        (e1-e5) + t1 t2 (1 + e1^2 - e4 - e1e5 + e2e6 - e6^2)`, with
        `e_k` in `(t1,t2,u1..u4)`.  It matches exact `phi_1` through
        degree 14.
      - No positive decomposition of the product class exists at
        degree 14 (8,412 `S_2 x S_4`-symmetric columns), nor at degree
        16 with monomials in `{0,2}^6` (15,489 columns).  This is the
        same obstruction as T1 at `s = 7`.
    - **Theorem R3-4 (FM-R3F, luna_max_mars; independently rerun by the
      main agent, `xmn/r3check.py`).**  FM3 at `r = 3` holds for every
      H-only word with at most four factors, for all degrees at once.
      - `A_3(kappa) = <h_kappa, Q_2^2>`, with
        `Q_3 = Q_2^2 = 64 + 64h_2 - 48h_1^2 + 16h_2^2 - 24h_2h_1^2 +
        9h_1^4`.  This extracts invariant Hilbert series on up to
        `s + 4` vertices.
      - Numerators:
        - `N^(3)_2 = 35 + 14m_2 - 21m_11 + m_4 - 4m_31 + 6m_22`;
        - `N^(3)_3 = N^(3)_2 + m_211`;
        - `N^(3)_4 = N^(3)_3 + 41m_1111 + 5m_3111 - 5m_2211`.
      - For `s = 4` the extraction uses the 8-vertex ring.  Mars
        determined the numerator through the a priori bound, degree 28.
        The main agent re-derived all three numerators from its own
        Weyl-character series (through degree 16, 14, 16).
      - Positive decompositions, all verified exactly by the main
        agent:
        - `s = 2`: 9 integer terms;
        - `s = 3`: 24 integer terms;
        - `s = 4`: 19 `S_4`-symmetrized terms with rational
          coefficients.
      - Boundary values `A_3 = 20, 8, 20` at `(2,2)`, `(2,1,1)`,
        `(1^4)`.  The forced zeros hold (`kappa_1 - |kappa|/2 >= 3`).
    - **Theorem S1R2-4 (FM-S1R2, luna_max_mars; main agent independent
      numerators and certificates, `xmn/s1r2_num.py`, `posdec_fix0.py`).**
      FM3 at `r = 2` holds for every word `hat S_p h_kappa` with at most
      four symmetric-power factors, for all `p` and all degrees.  This
      is the open-core sector "`r = 2` with one `hat S`", up to four
      factors.
      - Numerator: `N_s = (2+2t^2)(8T(0) + 4T(2) - 3T(1,1)) -
        t(8T(1) + 4T(1,2) - 3T(1,1,1))`, where `T(alpha)` extracts
        auxiliary degrees `alpha` from the invariant Hilbert series.
        - `N_1 = N_2 = B_s`, with
          `B_s = 10 + 2m_2 - 4m_11 - 7t m_1 + 9t^2 - t m_3 + t m_21 +
          3t^2 m_2 - t^2 m_11 - 3t^3 m_1 + t^4`.
        - `N_3 = B_3 + 9t m_111 - 2t^2 m_211 + 4t^3 m_111`.
        - `N_4 = N_3 + 6m_1111 - 4t m_2111 - t^2 m_1111 + t^2 m_3111 +
          t^2 m_2211 - t^3 m_2111 - t^4 m_1111`.
        - `N_4` comes from mars's exact extraction through the 8-vertex
          degree bound 28.
      - The main agent recomputed all four numerators from its own Weyl
        values of `<hat S_p h_kappa, Q_2>`, through degree 16, 16, 16
        and 14; they agree.
      - The main agent found its own `S_s`-symmetric exact certificates
        with 8, 14, 19 and 35 terms (`|m| <= 4`, all edge sets).  Each
        full expansion equals the numerator.  Mars's certificate rows
        for `s >= 2` were not printed, so the certificates are the
        main agent's.
    - **Level-1 small parts for fixed pairs (FM-SEC3, luna_max_saturn;
      conditional, accepted).**
      - For each fixed pair `(p,q)`, the words with `X_pq h_kappa` not
        genuine form a finite set.
        - Lemma BP bounds the parts by `p+q-2`.
        - Pure-power genuine seeds `k^(n_k)` bound each multiplicity:
          - `(2,2)`: `1^9`, `2^3`;
          - `(3,3)`: `1^17`, `2^8`, `3^4`, `4^3`;
          - `(4,4)`: `1^29`, `2^13`, `3^7`, `4^5`, `5^3`, `6^3`.
      - With exact character identities (FM53-cone certificates for
        `h_1`, `h_2^3`, `h_(1^3)`, `h_(1^4)`) and exhaustive finite
        checks, the level-1 sector is proved for ALL words for pairs
        `(2,2)`, `(3,2)`, `(4,2)`, `(3,3)`, `(4,3)`.
      - `(4,4)` leaves 34,411 residual words, not evaluated.
      - The all-ones seed threshold grows with the labels (9, 12, 15,
        17, 22, 29), so this does not give uniformity in `(p,q)`.
    - **Three labels times two symmetric powers at `r = 1` (FM-S3L1,
      luna_max_saturn; conditional, accepted).**  These are the
      Q3-relevant level-1 words with three plus labels `>= 2`.
      - Exact numerator in `X = (t1,t2,t3,u1,u2)`:
        `N = A1A2A3 - sum_i t_i prod_(j != i) A_j (e1 - e5) +
        sum_(i<j) t_i t_j A_k C_2 - t1 t2 t3 C_3`, where:
        - `A_i = 2 + 2t_i^2` and `C_2 = e1^2 + 1 - e4 - e1e5`;
        - `C_3 = 3m_1 + 5m_111 + 3m_21 + m_3 - 11m_11111 - 2m_2111 -
          m_22111 - m_31111` is the new 8-vertex term, taken through the
          Bruns--Herzog bound;
        - `N` has 98 terms and degree 10.
      - The main agent confirmed `N` from its own Weyl values through
        degree 12.
      - The degree-10 product class is infeasible, with an exact
        Farkas certificate using only rows of degree 0 and 2.  That
        certificate breaks for edge sets with more than 5 `t`-edges,
        so higher degree stays open.
      - **Theorem S3-2 (main agent; verifier
        `character_ring_iter/verify_shat_triple_two_factors.py`; ACCEPTED
        by FM-CHK24, luna_max_mercury: independent Weyl values on 6,188
        tuples, and exact re-expansion of the 70 terms).**
        FM3 at `r = 1` holds for every `hat S_p hat S_q hat S_w h_a h_b`,
        for all labels and degrees.
        - Degree 12 is still infeasible (phase-1 slack 0.138).
        - Degree 14 is feasible: 70 `S_3 x S_2`-symmetrized terms with
          nonnegative rational coefficients, reconstructed and verified
          exactly.
        - The 8-vertex term `C_3` rests on FM-S3L1's extraction through
          the degree bound.
    - **Degree of `N_7` (FM-GEN2, luna_max_venus; proved).**
      - Enumerate the Jozefiak--Pragacz--Weyman terms of `K_9`
        (Thms 3.3, 3.14; 16 terms, maximal `|lambda| = 36`).
      - This gives exactly
        `N_7 = 5 + s_2 - 3s_11 + 3s_1^4 - 5s_1^6 - s_21^4 + 5s_21^6 +
        s_2^21^4 - 3s_2^31^4 + 3s_2^51^2 - s_2^6 - 5s_2^7`, which is the
        numerical formula.  So `deg N_7 = 14`.
      - The positivity of `N_7 Pi_7` stays open.  No all-degree
        separator for the product class was found.
    - **Three factors at every level (FM-THREE, luna_max_jupiter;
      conditional).**
      - Explicit multiplicity formula (Pieri cells of
        `Sym^a (x) Sym^b` times the `Sym^c` weight pyramid, alternated
        over `W(C_2)`), matched on 9,007 slots.
      - Exact screen: `a+b+c <= 40`, `r <= 40`; 47,320 even-degree
        pairs.  There is no negative value, and zeros occur exactly in
        the forced range.  Odd total degree gives 0 identically.
      - Termwise positivity over Pieri cells fails, e.g. `(2,1,1)` at
        `r = 2` has cell contributions `20, -40, 60`.  The aggregate at
        `r >= 4` is open.  [Closed on 2026-09-30 by Theorem THREE (item (39)).]
    - **Theorem R3-5 (FM-R3F5, luna_max_mars; main agent reran the
      numerator and certificate, `xmn/r3s5check.py`).**  FM3 at `r = 3`,
      H-only, holds for every word with at most five factors.
      - `N^(3)_5 = N^(3)_4 - 24m_21111 - m_41111 - m_32111 +
        9m_22211 - 10m_22222`.
        - The 9-vertex term uses the JPW bound: degree of `N <= 40`.
        - Exact Weyl finite differences through degree 40 found no
          further terms.
        - The main agent matched it with its own series through degree
          14.
      - Positive certificate: 28 orbit-sum terms, verified exactly.
      - `s = 6`: the numerator is checked through degree 28 (support up
        to degree 14).  A restricted LP (`|m| <= 4`) is infeasible;
        this is open.
    - **Two-factor words at every level (FM-TWO, luna_max_jupiter;
      partial, accepted).**  Here
      `F(a,b;r) = 2r(2r+1) phi_r(h_a h_b)` is the alternating sum of the
      item-(35) weights over the Pieri support.  Put `s = b`,
      `d = (a-b)/2`, `N = 2r+1`, `M = r-1-d` and `k = max(0, M-s)`.
      - Exact reduction:
        `F = S_r(k,M) + N(N-1)(A_0 + (-1)^s A_1)`.
        - `S_r(k,M) = sum_(m=k..M) C(N,m+1) C(2r,m)
          (2(r-m)^2 + (2r-1)(r-m) - r)/(2r-1)` has positive terms.
        - `A_0, A_1 > 0` are explicit binomial products, and vanish when
          `k = 0`.
      - Consequences:
        - `F = 0` for `r <= d` (the forced zeros);
        - `F > 0` for `r > d` whenever `b` is even, or `k = 0`, or
          `A_1 <= A_0`.
      - Checks: the reduction against the original sum on 2,028 triples,
        and against direct Catalan moments.
      - **Odd `b`: PROVED (FM-TWO2, luna_max_jupiter; verifier rerun by
        the main agent).**
        - The sum telescopes:
          `S_r(k,M) = ((2r+1)/(2r-1)) (G(M+1) - G(k))`, with
          `G(m) = m(2r-m) C(2r,m)^2/(2r)`.
        - So `S_r - N(N-1)(A_1 - A_0)` is a positive factor times
          `P(x) = A x^2 - B x - C`, where `x = C(2r,k+s+1)/C(2r,k)`.
        - Jensen (convexity of `log((1+z)/(1-z))`) gives
          `x >= x_2 = 1 + 2(u+1)delta + (u+1)(2u+1)delta^2`, where
          `s = 2u+1` and `delta = 4(t+u)/(2k+2u+3)`, `t = r-k-s`.
        - `P(x_2)` and `P'(x_2)` are rational functions whose numerators
          have only positive coefficients (212 and 88 terms).  So
          `P(x) > 0`.
        - Main-agent independent check straight from the definition:
          `F(a,b;r)` over 13,230 cases (`a, b <= 40`, `r <= 30`).  There
          is no negative value, and `F = 0` exactly when `r <= d`.  `F`
          is also proportional to `<h_a h_b, Q_2^(r-1)>` on samples at
          `r = 2, 3, 4`.
    - **Theorem TWO (FM-TWO + FM-TWO2).**  FM3 holds for every two-factor
      H-only word `h_a h_b` at every level `r >= 1`.  It is `0` exactly
      on the support region `r <= (a-b)/2`, and positive otherwise.
    - **FM-EXP6 (venus): the Pieri-shift induction does not close.**
      - For `M = M' (x) V_k`,
        `Hom(U, M) = (+)_lambda E_lambda(M') (x) C_lambda(k)`,
        with `lambda` in `U (x) V_k`, and `M^G = E_(k,0)(M')`.
      - Old-slot maps become reduced slot maps of `M'` tensored with
        recoupling (6j) coefficients.  Last-slot maps are block diagonal
        in `lambda`.
      - So injectivity for `M` is equivalent to a shifted statement: one
        simultaneous injectivity over all Pieri cells `lambda` of `M'`.
        It does not follow from GSI for `M'`.
        - E.g. `M' = V_(2,0)`, `k = 2`: `Hom(U, M') = 0`, yet `M` has
          `m_U = 1`.
        - Several `lambda`-cells feed the same target, so cancellation
          after recoupling must be excluded jointly.
      - This is the linear analogue of the Pieri-induction obstruction of
        item (30).  Four or more slots need a different argument.
- **Rank-profile criterion (FM-INJ4, venus; proved).**
- **KILL (FM-GSI4, luna_max_venus): the rank-profile route at four
  factors.**
  - At `kappa = (6,6,6,4)`, `m_00 = 15`, `m_11 = 65`, `m_20 = 126`.
    The stacked map `f -> (A_1 f, ..., A_4 f)` into `(M^G)^4` has a
    kernel of dimension `>= 65 - 60 = 5`.
  - So the profile locus `alpha = 0`, `beta <= 4` has projective
    dimension `>= 4`, while the `d = 1` criterion needs `< 2`.  The
    bound fails.  The same happens at `(7,5,5,5)`.
  - T1 there is 6, and GSI itself is not refuted.  The generic map is
    `195 -> 201` with slack 6.
  - Main-agent direct test (`inj/gsi_6664.log`, mod 1,000,003):
    `(6,6,6,4)` has `rank(T) = 195/195`, INJECTIVE, in 1,509 s, and
    `(7,5,5,5)` has `195/195`, INJECTIVE, in 1,217 s.
    - The `B`-channel absorbs the 5-dimensional kernel of the stacked
      `A`-map.  So the obstruction is to the rank-profile criterion only.
    - A proof of GSI has to use `A` and `B` jointly: `alpha(f) = 0` loci
      of positive dimension really occur.
  - FM-GSI5 (luna_max_jupiter): a joint criterion; conditional.
    - Exact `K_A = ker(A_1..A_4)` on all 973 four-part words of size
      `<= 26`: `dim K_A` is 0 for 825 words, 1 for 76, 3 for 42, 6 for
      20, 10 for 8, and 15 for 2.
      - At `(6,6,6,4)`, `dim K_A = 6`; at `(2,2,2,2)` it is 1.
    - Sufficient criterion (proved):
      1. `B` is injective on `K_A^3`;
      2. `B` is injective on the residual `A`-kernel `N` modulo
         `B(K_A^3)`.
    - The simple split (`N = 0`) fails by dimension on 338 words, first
      at `(2,2,1,1)`.  At `(6,6,6,4)`, `dim N >= 102`.
    - So the joint condition is a reformulation, not a proof.
  - Boundary words `(1^4)`, `(2,2,1,1)`, `(3,1,1,1)`, `(2^4)`,
    `(4,2,1,1)` have no such forced kernel.
- **Covariant model of the slot maps (FM-DEG1, luna_max_jupiter;
  partial).**
  - Covariants: `Hom(U, M) = (U (x) C[z_i])^Sp4`.  This is generated over
    `R = C[omega_ij]/(Pf_6)` by `C_ij = X_(z_i, z_j) - omega_ij I / 2`.
  - `A_i` and `B_i` are first-order operators given by the product rules
    `D_X^(i) omega_rs = [i=r] omega(X z_r, z_s) + [i=s] omega(z_r, X z_s)`
    and the analogue for `C_jk`.  Explicitly,
    `A_i(c C_jk) = (1/2) sum_(l != i) (d c / d omega_(il)) tr(C_(il) C_jk)
    + (5/2)([i=j] + [i=k]) c omega_jk`, using
    `tr(C_ij C_kl) = 2(omega_ik omega_jl - omega_il omega_jk
    - omega_ij omega_kl / 2)`.
    - CORRECTED endpoint coefficient: `5/2`, not the `1/2` originally
      printed by FM-DEG1.
      - FM-EXP4 (venus) contracted exactly in `slotHWg.py`'s conventions
        and got `A_0(F_12) = 3 omega^2` at `(2,2)`.
      - Independently, with `1/2` all three `A_i` would be proportional on
        `H` at `(2,1,1)`, giving joint-`A` rank 1.  The measured joint-`A`
        rank there (`hallstats.py`) is `2 = m_U`.
      - The `1/2` came from missing the `C`-factor contraction, which
        contributes `5/2` rather than `1/2`.
  - The naive covariant basis (3-nesting-free graphs with one marked edge)
    overcounts.  It gives 6 vs `m_U = 5` at `(1^4)` (42 vs 30 at `(1^6)`),
    because of the syzygy
    `C12 w34 + C34 w12 - C13 w24 - C24 w13 + C14 w23 + C23 w14 = 0`.
  - A module term order that makes `C14 w23` leading fixes `(1^4)`
    locally.
  - **FM-DEG2 (jupiter): module bases and graded maps.**
    - Presentations: `(Sym^2 W (x) S)^G = R{q_ij = z_i z_j : i <= j}`.
      - `U` relations are the Fierz syzygies
        `w_kl C_ij + w_ij C_kl - w_jl C_ik - w_ik C_jl + w_jk C_il
        + w_il C_jk = 0`.
      - `Sym^2 W` relations are the five-vector relations
        `sum_r (-1)^(r-1) Pf(w_(J \ j_r)) q_(j_r k) = 0`.
    - Use the anti-diagonal order on edges (3-nesting-free invariant
      basis) and position-over-term module orders.  The standard-monomial
      counts then equal `m_U = 1, 2, 5, 3, 3, 7, 30` and
      `m_ad = 1, 3, 6, 6, 6, 10, 40` on `(2,2), (2,1,1), (1^4), (2,2,2),
      (3,1,1,1), (2,2,1,1), (1^6)`.  The relation ranks account for the
      full kernels there.
    - Graph moves: `gr A_i` consumes an edge `(i,l)` and the marked edge
      `(j,k)`, and produces the reconnections `(i,j)(l,k)`, `(i,k)(l,j)`,
      `(i,l)(j,k)` with coefficients `+1, -1, -1/2`.  There is also an
      endpoint term, `(5/2) c w_jk` after the correction above (printed
      `1/2` in FM-DEG2).  `gr B_i` gives a marked symmetric pair `q_ab`.
    - The generic `3 -> 5 + 1` combination of the graded maps has full
      column rank on all seven words.
    - Open:
      - a common flat (Rees) filtration showing that these leading terms
        are the associated graded of `A_i`, `B_i`.  Semicontinuity then
        gives injectivity of `T`.
      - a uniform argument.  Hall on the support graph is not enough,
        since coefficients are shared across rows; a non-cancelling
        matching is needed.
  - **FM-DEG3 (jupiter): obstruction to the naive Rees degeneration.**
    - The slot-multidegree weight `w(omega_ij) = w(C_ij) = w(q_ij) =
      e_i + e_j` preserves the full reconnection stencil, and `A_i`, `B_i`
      are filtered for it.  But it has one level within each `kappa`, so
      it degenerates nothing.
    - Any scalar edge weight that ties the three pairings of every four
      vertices has the form `w_ij = a_i + a_j`.  Such weights give all
      Pfaffian terms equal weight and cannot select the anti-diagonal
      (3-nesting) initial term.
    - So a degeneration that realizes the standard-monomial basis must use
      non-additive weights.  Then `A_i` is filtered only up to a shift
      `delta`, and `gr A_i` keeps only the top-weight reconnection.
    - Open: whether these partial graded maps stay injective.
    - **KILL (FM-DEG4, jupiter): uniform-shift Rees degeneration.**
      - Setup: weight `w(omega_ij) = (j-i)^2`, which makes the nested
        matching the unique top term of each Pfaffian.  Module positions
        are ranked by span, with strict position priority.
      - `(2,2)` passes (3/3).
      - At `(2,1,1)` the second source generator `omega_12 C_13` has every
        image below the top gain.  So `gr A_i` and `gr B_i` kill it, and
        the graded rank is 3/6.
      - The original map is injective there (6/6), so this kills only the
        degeneration route with this order and uniform shifts.  A working
        degeneration would need source-dependent shifts, which is not a
        single Rees family.
    - Precise graph form of GSI at `r = 2`: generic combinations of the
      `A_i` edge-switch and `B_i` marked-pair moves, from three copies of
      the Fierz-reduced marked 3-nesting-free basis into five copies of
      the invariant basis plus one symmetric-pair basis, are injective.
- **Literature (FM-LIT4, luna_max_neptune): no known theorem closes T1.**
  - Confirmed exactly: `8H(kappa) + 4H(kappa,2) - 3H(kappa,1,1) =
    5m_00 - 3m_11 + m_20`, where `H` is the multigraded Hilbert function
    of the rank-`<= 4` skew Pfaffian ring.
    - Checked on `(1^4)`, `(2,2,2)`, `(3,1,1,1)`, `(1^6)`: margins
      6, 2, 2, 20.  The free ring gives -15 at `(1^6)`.
    - Gorenstein structure alone cannot give the inequality.
    - **FM-PF1 (luna_max_jupiter): multigraph form, one class killed.**
      - Put the two new vertices last.  The edge `{s+1, s+2}` can
        neither nest nor be nested, so `H(kappa,1,1) = H(kappa) +
        P11(kappa)`.  Here `P11` counts 3-nesting-free graphs with
        pendant edges at `s+1` and `s+2`, and `H(kappa,2) = P2(kappa)`.
        So T1 is `3 P11 <= 5 H + 4 P2`.
      - Exact enumeration on all 139 partitions with `|kappa| <= 10`,
        and on all 256 part orders of size `<= 8`, gives `H = m_(0,0)`,
        `P11 = m_(1,1) + m_(2,0)` and `P2 = m_(2,0)`, independent of the
        order.
      - KILL (direct merge/slide maps).  At `(1^6)`, the two `P11`
        graphs with old edges `(1,6), (2,5)` and pendant edges
        `(3,7), (4,8)` or `(4,7), (3,8)` both slide to the 3-nesting
        `(1,6), (2,5), (3,4)`.  Both also merge to the same `P2` graph.
        So six source objects meet four target copies.
      - A class with rerouting moves was not specified precisely enough
        to test.
      - KILL (FM-PF2, luna_max_jupiter): merge/slide with canonical
        local repair.
        - The repair swaps right endpoints of the innermost nested pair,
          shortest repair path, and any route choice or copy label is
          allowed.
        - A Hall obstruction appears at `(3,1,1,1)`: 18 tagged sources,
          all sliding to the same star, can reach only 17 targets.  The
          maximum matching is 26 of 27.
        - All words of total degree `<= 5` match.
  - Enright--Willenbring (Ann. Math. 159 (2004)):
    - Thm 3(i) gives Howe duality `(Sp4, so*(2s))` on `P(M_(4 x s))` in
      the full range.  So FM3 H-only at level `r` is weight-positivity of
      `sum_mu c_mu(r) L_mu`.
    - Their Thm 2 (generalized Verma resolutions of each `L_mu`) and Thm 5
      (Hilbert polynomials) do not control the signed combination.
  - Jonsson--Welker (math/0601335):
    - Prop 3.1 is the 3-nesting-free Herzog--Trung basis.
    - Thm 2.1 is a second term order (circular distance plus revlex), with
      3-crossing-free long edges.  It has a different symmetry; relevant
      to FM-DEG1.
  - Ghorpade--Krattenthaler (Thm 4.1, Cor 4.3) and De Negri give ordinary
    Hilbert series only.
- **Kostka literature for MP_2 (FM-LIT5, neptune): no applicable theorem;
  one local injection killed.**
  - None of these supplies the signed shape sum:
    - White / Fayers Prop 1.2 (monotonicity in content, fixed shape);
    - Bender--Knuth and RSK (Knuth Thms 2, 7);
    - LPP Thms 4--5 (Schur-positive outputs; `F_2` is not
      Schur-positive: `F_(2,4) = 5s_(1^4) + 2s_(22) + s_(31) - 2s_(211)`);
    - Lam--Pylyavskyy cell transfer Thm 3.6 (acts on pairs of tableaux);
    - Johnston et al. Prop 3.5 (a different box move);
    - Kirillov Thm 3.3.
  - Exact Kostka screen: MP_2 holds for all 4,095 compositions of size
    `<= 12`, with no negative value.
  - **KILL (same-column-count local injection).**  Map `(k+1,k,1)+j^4` to
    `2 x (k+1,k+1)+j^4`, `(k+2,k)+j^4` and `(k,k,2)+j^4`.
    - This fails from degree 6.  At `(2,1^4)`, `k = 2`:
      `2 K_(321) = 16 > 2 K_(33) + K_(42) + K_(222) = 6 + 6 + 2 = 14`.
    - The full margin there is still 8.  The extra capacity comes from
      shapes with one more full column, `(2,2,1,1)` and `(3,1,1,1)`.
    - So a tableau injection must also move boxes into row 4, changing the
      column count.
- **Rank-profile criterion (FM-INJ4, venus; proved).**
  - Generic `T` is injective if every profile stratum
    `G_(d, alpha, beta)`, `1 <= d <= min(3, m_U)`, satisfies
    `dim G + 3d - 1 < 5 alpha + beta`.
    - Here `G` is the set of subspaces `H'` of `Hom(U,M)` with
      `dim H' = d`, `alpha = dim sum_i A_i(H')` and
      `beta = dim sum_i B_i(H')`.
    - The proof is an incidence-dimension count.
  - Checked by hand at `(2,2,2)` and `(3,1,1,1)`, from exact minors: joint
    `A` injective, individual `B_i` of rank 3, pairwise `B` spans of 5.
    So T1 holds there by this route as well.
  - The global profile bound is open; it is equivalent in strength to
    GSI at `r = 2`.
- **Independent check of the GSI screens (FM-GSIX, luna_max_mars,
  independent code).**
  - Reproduced exactly, with no shared code:
    - `r = 2`: `(2,2,2)` 9, `(3,1,1,1)` 9, `(1^4)` 15;
    - `r = 3`: `(2,1,1)` 80, `(3,1,1,1)` 145, `(2^4)` 690;
    - `r = 4`: `(2,1,1)` 966.
  - New injective cases:
    - `r = 3`: `(4,4,1,1)` 320/320, `(5,3,1,1)`, `(4,3,3)`, `(4,4,2)`,
      `(5,4,1)`, `(5,3,2)`, `(5,5)`;
    - `r = 5`: `(5,1)` 231/231.
  - Main agent, the tightest cases found by mars's multiplicity screen:
    - `r = 3`, `(2^4,1,1)`: 6275/6275, slack 40 (0.64%);
    - `r = 3`, `(3,2,2,1,1,1)`: 4155/4155, slack 28 (0.67%);
    - `r = 5`: `(3,1)` 1650/1650, `(4,2)` 1881/1881, `(4,1,1)` 3762/3762,
      `(2,2)` 5808/5808, `(3,3)` 6963/6963, `(2,1,1)` 11616/11616 and
      `(3,2,1)` 16302/16302, all injective.  The smallest relative slack
      is 3.1% at `(2,1,1)`.
  - A single-slot cyclic rule in edge order at `r = 3` is NOT injective at
    `(2^4)`: 687/690 on two primes.  Single-slot patterns do not carry over
    from `r = 2` automatically.
- **Structure of the `r = 2` slot maps (`hallstats.py`, 17 words to size
  12).**
  - The joint map `(B_1, ..., B_s): Hom(U,M) -> Hom(ad,M)^s` is injective
    in every case.
  - The joint `A` map has a kernel only at `(2^4)`, `(3,3,2,2)`,
    `(4,4,2,2)` (dimension 1) and `(3^4)` (dimension 3).
  - `sum_i A_i(H) = M^Sp4` always.  `sum_i B_i(H) = Hom(ad, M)` except at
    `(2,2,2)` (5 of 6).
  - On `ker(joint A)`, the `B`-images span `4 dim ker`.  This is the linear
    Hall condition there, with slack `dim ker`.
- **Method note (lean certificate).**
  - The dense stage is reduced to the kernel of one simple-root raising
    operator, computed block by block:
    - `e_(a1) = x0 d/dx1 - x3 d/dx2` preserves each slot's
      `(e0+e1, e2+e3)`;
    - `e_(a2) = x1 d/dx3` preserves `(e0, e2)`.
  - Blocks are tensor products of `sl_2` modules.  The kernel has
    dimension `dim M_mu - dim M_(mu+alpha)`, about `1/9` of `M_(1,1)` for
    `alpha = alpha_1`.
  - The multiplicities come from the Weyl character and are cross-checked
    against the nullity.
  - Images are compressed by sparse random projections.  These only lower
    the rank, so full rank still certifies.
  - The negative control reproduces the fixed-power failure at `(2,2,2)`.
- **What this gives the full cone.**  For the H-only sector, GSI is one
  uniform statement:
  - it covers every level;
  - it holds in every case tested (`r = 2, 3, 4`);
  - it implies FM3 there.
  For `hat S` words the lift works with one label and fails, as built,
  with two (see the kill above).  So GSI does not yet cover the full cone.
  The open step is a uniform proof of generic injectivity, and it must be
  constructive.  Lemma JB (joint `B` injective) is a first structural
  piece; FM-INJ4 (luna_max_venus) is assigned to it.
- **Hall/rank analysis (FM-RADO, luna_max_saturn; main-agent check).**
  - Proved: injectivity implies the linear Hall condition (LH).
  - Proved: LH holds iff `T`'s block space has full non-commutative rank.
    - Fortin--Reutenauer Thm 1 gives the shrunk-subspace formula.
    - The deficiency `dim Z - dim A(Z)` is supermodular and invariant
      under `prod GL(c_mu)`.  So a maximal shrunk subspace can be taken
      to be `(+) X_mu (x) C^(c_mu)`.
  - LH does not give commutative generic injectivity in general.  The
    `3 x 3` skew-symmetric space at copy number 1 satisfies LH and has
    rank 2.
  - Derksen--Makam Thm 1.8 closes the gap only for uniform copy scaling
    `c = d c^0` with `d >= m_0 - 1`.  That does not cover the T1 profile
    `(3; 5, 1)`.
  - **Consequence (main agent).**  FM3 is LH at `X_mu = H_mu`, and
    Conjecture T is LH on `{0, H_mu}`.  So "LH implies injectivity" is no
    route to FM3.  GSI helps only through a constructive injectivity proof,
    for example:
    - a coefficient specialization that is block-triangular for a
      filtration of `Hom(V_mu, M)` (by slot support or by degree), with
      injective graded pieces;
    - an induction on the last tensor factor.
    The screens show GSI is true, and LH with it, which rules out a
    subspace obstruction in these cases.  They do not supply the proof.

**Ledger for this section (2026-09-28).**

- **Proved (checked by FM-CHK18).**
  - FM53, reproved via Kostant and extended to every `GL(4)`-module `E`
    (item (5)).
  - The `y -> -y` duality of words (item (8)).
  - qFM3_2 in these cases (item (9)): one part; two parts; the support
    region.
  - The all-ones product formula (item (1)).
  - Item (35): the closed form of the FM12 weights (GKS reflection), the
    real-level product formula for irreducible characters, and the
    forced-zero divisibility for every FM3 word (FM-P35, reviewed).
  - Item (38), each checked by FM-CHK20:
    - Lemma JB (each `B_i` injective), hence `m_(1,1) <= m_(2,0)` for
      every nonempty product;
    - LH iff full non-commutative rank for copy-blown-up block spaces
      (FM-RADO);
    - the dimension-0 constraint, which rules out whole-module proofs of
      T1.
  - Item (38), added 2026-09-29:
    - Theorems T1-5 and T1-6: T1 (FM3, `r = 2`, H-only) for every word
      with at most six factors (T1-6 ACCEPTED by FM-CHK-T16).
    - Theorems R3-4 and R3-5: FM3 at `r = 3`, H-only, for every word
      with at most five factors (mars; main-agent rerun).
    - T1 for seven-factor words whose smallest part is `<= 4` (main
      agent, slicing, exact certificates).
    - Theorem TWO: FM3 for every two-factor H-only word at every level
      (jupiter FM-TWO, FM-TWO2; main-agent rerun and direct check).
    - Theorem S1R2-4: FM3 at `r = 2` for one `hat S` label times at most
      four symmetric powers, all labels (mars; main-agent independent
      numerators and certificates).
    - Theorem S3-2: FM3 at `r = 1` for three `hat S` labels times two
      symmetric powers, all labels (Q3-relevant; main agent, exact
      verifier).
    - Level-1 two-label sector complete for the pairs `(2,2)`, `(3,2)`,
      `(4,2)`, `(3,3)`, `(4,3)` (saturn FM-SEC3; exact identities and
      finite checks).
    - Theorems S2-2 and S2-3: FM3 at `r = 1` for `hat S_p hat S_q`
      times at most three symmetric powers, for all labels (main agent;
      exact verifiers; ACCEPTED by FM-CHK-S2).  The proof uses generating functions from the
      `Sp(4)` fundamental theorems with the Buchsbaum--Eisenbud resolution,
      plus an exact positive decomposition.  The verifiers are in
      `character_ring_iter/`.
    - Theorem (three slots, FM-EXP5): `T_exp` is injective for every
      three-factor word, with a verified certificate.
    - The level-1 identity `hat S_p hat S_q = sum_m hat S_m + X_pq`,
      with `X_pq` a diamond of four irreducibles (FM-CHK22 accepted).
    - Lemma BP and its corollary (main agent; ACCEPTED by FM-CHK-BP.  Mars
      also shows the threshold cannot drop to `sum p_j - 2`: labels
      `(2,2,2)` at `k = 4` give `-3, -3, -1` at `V_11, V_33, V_55`).  At
      `r = 1`, FM3 holds for `prod_j hat S_(p_j) h_kappa` whenever some
      part of `kappa` is `>= sum p_j - 1`.  This holds for any number of
      labels.
    - Labels `(2,2)` at `r = 1`: every `kappa`.  The ingredients are
      Lemma BP, two exact FM53-cone identities (for `h_1` and `h_2^3`),
      and direct values at `()`, `(2)`, `(2,2)`.
- **Proved strata of FM3 (item (38) generating-function method and
  related results; state at 2026-09-29 13:30).**

  | level | sector | proved for |
  |---|---|---|
  | `r = 1` | H-only | all words (invariant count) |
  | `r = 1` | one `hat S` | all words (FM53) |
  | `r = 1` | any labels `p_1..p_l` | words with a part `>= sum p_j - 1` (Lemma BP) |
  | `r = 1` | two labels | `<= 3` symmetric powers, all labels (S2-2, S2-3); all words for pairs `(2,2), (3,2), (4,2), (3,3), (4,3)` (FM-SEC3) |
  | `r = 1` | three labels | two symmetric powers, all labels (S3-2) |
  | `r = 2` | H-only (T1) | `<= 6` factors (T1-5, T1-6); 7 factors with smallest part `<= 4` (slicing); every `d = 1` word (Lemma D1) |
  | `r = 2` | one `hat S` | `<= 4` symmetric powers, all labels (S1R2-4) |
  | `r = 3` | H-only | `<= 5` factors (R3-4, R3-5); `h_s h_t h_1^a`, all `s, t, a` (R3-TWO, item (39)) |
  | `r = 3` | one `hat S` | `hat S_p h_1^a`, all `p, a` (R3-S1, item (39)) |
  | all `r` | H-only | two factors (Theorem TWO) |
  | all `r` | one label | `hat S_q h_1^a` and `h_s h_1^a`, all `q, s, a` (Theorem OL, item (39)) |
  | all `r` | two labels | all three sign patterns when the suffix `a` is `<= 2` or within 1 of the weight exponent `e` (item (39), strip of (E)) |
  | all `r` | H-only | three factors, no suffix: `h_u h_v h_w` (Theorem THREE); three factors, any suffix, `min label >= a+2r-4` or `u-v+w >= a+2r-3` (Theorem LL); any number of factors, all parts even, no suffix (item (39)) |
  | all `r` | H-only | three factors with `|a - (2r-3)| <= 1`, all labels (Theorems DS, AS); T3R regions (item (39)) |
  | all `r` | H-only | three factors in the G0 branch `gamma >= a+2r-2` under the one-sign-block, root-crossing (G3X) or odd-`C` (G3O) conditions, and parts of `gamma = N` (Theorems G3, G3X, G3O, item (39)) |
  | all `r` | H-only | three factors, whole G0 branch when `(a-2r+2)^2 <= 8r-9` (Theorem G0B); whole G0 branch conditional on (E) (Theorem G0E) |
  | all `r` | H-only | four factors on the outer-pairing support region, conditional on (E) (Theorem G0E4; strictly extends LL4) |
  | `r <= 8` | H-only | three factors, both branches, `a <= 40`, `w <= 14`, `u <= v+w+1` (finite box; every word has a positive definite binary form) |
  | all `r` | H-only | three factors with equal largest labels in G0: `h_u^2 h_w h_1^a` for `2u >= a+2r+w-2` (centered-window theorem, FM-MECH25); the whole G0 branch for every `r <= 100`, `a <= 150` (every one of 130.9M windows covered by proved criteria) |
  | all `r` | two labels | (E) at `q = 1` (gap `i = j+2`) for all `a, e` (Theorem EQ1); `q = 2` with the minus-W sign (OL on `(1-z^2)P`) |
  | all `r` | two labels (all three sign patterns) | every word whose kernel row has `a, e <= 80` (each of 6.3M (E) pairs covered by a proved criterion: binary forms, LD/metric chords, RF, EQ1, strips); sampled to `a, e <= 300` with no gap |
  | all `r` | H-only | any number `m <= 2r` of factors with spread labels: `lambda(S)+lambda(S^c) > a+2r-m` for every split (Theorem LLm; includes LL, LL4) |
  | all `r` | every consumer word (h and hat S) | spread labels over every two-sided assignment (Theorem LLm-S) |
  | all `r` | every word | suffix `a >= (2r+3) Lambda(w) - 2r - 4` (Theorem LS); `2r + a + t + 6 >= 7L(w)`, `L` an additive quartic in the labels (Theorem LR4) |
  | all `r` | every word with labels `<= 2` | all (FM-MECH41, positive Walsh realization, B) — added 2026-10-01 |
  | all `r` | labels `<= 2` plus one label `n <= 6` | all (FM-MECH44, FM-MECH45) — added 2026-10-01 |
  | all `r` | every word with labels `<= 3` | all (FM-MECH47, exponential ray suppression plus one fixed box; FM-CHK56) — added 2026-10-01 |
  | all `r` | every word with labels `<= 4` | all (FM-MECH48, uniform contraction lemma plus one box of 54,387,664 values; FM-CHK57) — added 2026-10-01 |
  | all `r` | every word with labels `<= 5` | all (FM-SEC135 suppression for `H >= 68` plus the exact box of 227,336,512,420 values, FM-SEC139; FM-CHK79 ACCEPT) — added 2026-10-01 |
  | all `r` | every word, all labels | weighted count of factors of label `>= 3` at least `T_0(k) = O(log k)`, `k` the largest label; any number of 1's and 2's (FM-MECH49; FM-CHK59); constant weighted cutoff `T >= 2^21` (ADV-2; FM-CHK61) — added 2026-10-01 |
  | all `r` | two labels, one of them 2 | all labels, both signs: (E) at `q = 2` (FM-MECH50, Theorem Q2+; FM-CHK58) — added 2026-10-01 |
  | all `r` | EVERY list at list distance `<= 4` | all (FM-MECH38 for `<= 2`; FM-MECH68 for 3, 4; FM-CHK67) — added 2026-10-01 |
  | all `r` | EVERY list at list distance 5 or 6 | all (FM-MECH71; FM-CHK69 item 1 accepted, items 2 and 3 repaired or withdrawn in FM-MECH80, check FM-CHK74 running); every distance `delta` once `max(t,b,k) >= (100 delta)^(300 delta^2)` (FM-MECH71) — added 2026-10-01 |
  | all `r` | every list, any distance `delta` | `N >= 432(delta+1)^2` (FM-MECH85; FM-CHK76 ACCEPT), or `k >= 40 delta` labels `>= 3` (FM-MECH85), or one non-distinguished label above `delta` with any smaller background (FM-MECH85); earlier `N >= 3200(delta+1)^3 + 760(delta+1)^2`, where `N = t + 2b + k` counts the non-distinguished factors after splitting non-distinguished `-2`; also the integer test (1), and `N >= 380(delta+1)^2` when the signed label-3 count is 0 (FM-MECH77; FM-CHK70 ACCEPT) — added 2026-10-01 |
| all `r` | two core labels on a {1,2} background, `b = 1` copy of `hat S_2` | balanced backgrounds `a = e`, both signs, all labels (FM-MECH76 Thm 2; FM-CHK73); every `b`: the central band `4|Delta| <= sigma`, `sigma/4 + 2 <= n - m - 2b <= sigma/2`, `m >= 35 + 12b` (FM-MECH76 Thm 3) — added 2026-10-01 |
| all `r` | two labels on the fundamental background ((E), `q >= 3`, both signs) | every row with anchor `j >= N/2` (the consumer's range): the regions of FM-MECH55..61, 66, 70, 73, 74 and the last anchors `X > 7 omega/10` (FM-MECH79); FM-CHK71 ACCEPT with that scope — added 2026-10-01 |
| all `r` | one label `p <= 10` on any {1,2} background | all ratios and signs (FM-MECH78: Bernstein cutoffs plus 3,510,309 exact representatives; check FM-CHK72 running) — added 2026-10-01 |
| all `r` | EVERY list at list distance 7 | all (FM-MECH80: exact distance-7 polynomial, tails `t >= 249` or `b >= 137` or `k >= 333`, Schur reduction, all-allocation finite part; FM-CHK74 items 1-3 and FM-CHK74b ACCEPT) — added 2026-10-01 |
| all `r` | every word of length `L >= 4` | minimum label `>= 2^(L-2) - 3`, and more generally `beta <= 1` (FM-MECH83) — added 2026-10-01 |
| all `r` | `h >= 2` non-distinguished labels above the distance `delta`, any smaller labels | background weight `D <= min(delta + h - 2, 2 delta + 2)` (FM-MECH83) — added 2026-10-01 |
| all `r` | one label `p` on a balanced or adjacent background (`|a - e| <= 1`) with `b` copies of `hat S_2` | `b <= min(a,e)` and `s^2 >= 4b(min(a,e) + 2)`, `p = 2s + |a - e|` (FM-MECH82) — added 2026-10-01 |
| all `r` | at most two labels `>= 5` (any size, any sign) on a labels-`<= 4` background of weight `D <= 60`; at most three for `D <= 48`; any number of saturated labels with `u <= 0, 1, 2` unsaturated for `D <= 60, 56, 48` | all (FM-MECH87: 5,285,147 + 1,388,843 profiles, 4.3e9 exact three-label kernels) — added 2026-10-01 |
| all `r` | labels `<= 4` plus saturated labels, `h >= 2` | `delta >= ceil(15D/16)`, or the count test `16^(D-delta-1) R_0/lambda <= 1` (FM-MECH87) — added 2026-10-01 |
| all `r` | every list with exactly two odd labels, any even background | reduced to lists with fewer factors (two-odd fusion, FM-MECH143); a minimal counterexample is pair-free with 0 or `>= 4` odd labels — added 2026-10-02 |
| `r = 2` | `-a,-b` odd, `-c,-d` even, `+2^t`, even `+e_i >= 4` | `d >= a + b + sum e_i`, every `t` (FM-MECH143) — added 2026-10-02 |
| all `r` | four largest labels `q + t_i` (pair-free quartet) over a word of weight `D` | `q >= D + M + 1`, `M` the quartet spread: reduced to fewer factors (FM-MECH142); all 7- and 8-factor such words; every all-minus word of 8 consecutive labels — added 2026-10-02 |
| all `r` | EVERY list with at most seven factors | all labels, all signs (FM-MECH153 Theorems 2, 3 with Prop. 1C and Corollaries 23A9ZZ10, 5A7B51) — added 2026-10-02 |
| all `r` | eight factors with no negative 3|5 term (e.g. every all-odd eight-factor list) | all (FM-MECH153 Lemma 4, Corollary 5) — added 2026-10-02 |
| all `r` | eight factors: all labels `>= 3`, or all `>= 2` with an odd label | all (FM-MECH155 with FM-MECH153 and two-odd fusion) — added 2026-10-02 |
| all `r` | nine factors with minimum label `>= 28` | all (FM-MECH159) — added 2026-10-02 |
| all `r` | EVERY list with at most ELEVEN factors | all labels, all signs (FM-MECH173 with FM-MECH169; independently checked: FM-CHK109, FM-CHK109b; finite box rerun by main agent) — added 2026-10-02 |
| all `r` | EVERY list with at most TEN factors | all labels, all signs (FM-MECH169 with FM-MECH164; FM-CHK108 ACCEPT) — added 2026-10-02 |
| all `r` | EVERY list with at most NINE factors | all labels, all signs (FM-MECH164 with FM-MECH160; FM-CHK106 ACCEPT) — added 2026-10-02 |
| all `r` | ten factors with minimum label `>= 16` | all (FM-MECH164) — added 2026-10-02 |
| all `r` | EVERY list with at most EIGHT factors | all labels, all signs (FM-MECH160 with FM-MECH153, FM-MECH155, two-odd fusion; FM-CHK99 and FM-CHK102 ACCEPT) — added 2026-10-02 |
| all `r` | `L` = 9..16 factors with minimum label at least 8, 10, 13, 17, 22, 28, 35, 45 respectively | all larger labels, all signs (FM-MECH166, direct positivity; FM-CHK107 ACCEPT) — added 2026-10-02 |
| all `r` | any large labels on a labels-`<= 4` background | `B >= max(384 sum (n_i+1)^2, 2^17)` (FM-MECH63/69; FM-CHK65, FM-CHK68) — added 2026-10-01 |
  | all `r` | labels `<= 2` plus one arbitrary label | `b <= 2` copies of `hat S_2` (OL, Q2+, B2), and every `b` once `nu + b >= 16(p+1)^4 + 2` (FM-MECH51; FM-CHK60, `b = 0` via OL), improved to `nu + b >= 4p^2 - 2` (FM-MECH52); `min(a,e) <= 1` for every `b` (FM-MECH52); distance 3 for every `b` (ADV-1); distances 4, 5 for every `b` (FM-MECH53) — added 2026-10-01 |

  - **Checkpoint 2026-10-01 21:30 (main agent).**  FM3 full cone: NOT
    finished.  Proved (independent checks in brackets):
    - every list at distance `<= 7` (FM-MECH68 [CHK67], FM-MECH71
      [CHK69, repaired in FM-MECH80], FM-MECH80 [CHK74 items 1-3; 580-row
      check, FM-CHK74b ACCEPT at 22:00]);
    - every list with `N >= 432(delta+1)^2` or `k >= 40 delta` (FM-MECH85
      [CHK76]);
    - every list with a non-distinguished label above `delta` (`h >= 1`;
      FM-MECH85 [CHK76], FM-MECH93);
    - every word with labels `<= 5` (FM-SEC135 + FM-SEC139 [CHK79]);
    - (E) on the fundamental background (FM-MECH79 [CHK71]); H = 1 for labels
      `<= 10` (FM-MECH78 [CHK72]); bounded small backgrounds (FM-MECH87
      [CHK78 pending]); many strata of H = 1 and two cores (FM-MECH76/81/82/
      84/86/89 [CHK73, CHK75, CHK77 pending]).
    - Remaining cone, exactly: the unsaturated region `h = 0`.  All
      non-distinguished labels are `<= delta`, with `delta >= 8`,
      `N < 432(delta+1)^2` and `k < 40 delta`.  The consumer is
      `p_delta >= p_(delta-1)` for the coefficients of
      `P(q) = E_eta prod_i (G_(n_i)(q) + eps_i q^(n_i/2) U_(n_i)(eta))`.
      By `k`:
      - `k = 0`: the H = 1 layer (FM-MECH92);
      - `k = 1`: the two-core layer (FM-MECH91, FM-MECH94);
      - `k >= 2` (FM-MECH95).
    - Update 22:40 (FM-MECH102): by the pair reduction it suffices to
      prove FM3 for PAIR-FREE lists (one sign per label value) with
      `k >= 2` non-distinguished labels `>= 3`.  Pair-free lists with
      `k <= 1` are covered by OL, FM-MECH52, (E) and FM-MECH64 (FM-MECH64:
      FM-CHK83 ACCEPT; the reduction itself: FM-CHK82 pending).
    - Update 23:10 (FM-MECH106): reduction to two monotonicity lemmas (M+)
      and (M-).  CORRECTED 23:45: both are false as stated (see the
      correction entry; part of the screening was vacuous by parity).  Only
      (M+even) and the T2 case survive.  Route B (pair-free `k >= 2`
      directly) is the main route.
    - Update 00:50 (FM-SEC148/151): a single conjectural lemma, LEMMA R
      (the Rule 1' removal never increases `g_p`), would give the whole cone
      with FM-MECH102.  KILLED 01:15 (counterexamples at W = 128, 131; see
      the KILL entry); only the existence form (D) survives.
    - Update 01:45 (FM-SEC156): new candidate LEMMA W (Rule W, by class
      weight), exhaustive through W = 40, structured families to W = 2030,
      adversarial to W = 120, no failure.  KILLED 02:40 (W = 427 and two
      more); Rule M and uniform averaging also killed.  Weighted averages of
      (D) killed later.
    - Update 03:45 (FM-SEC160): conjecture (S) killed; LEMMA WM (one of the
      Rule W or Rule M removals is monotone on the residual) survives with
      margin >= 0.648 adversarially to W = 700.  KILLED 04:25 at W = 1283
      (FM-MECH137).  Only the existence statement (D) survives.  It survives an exhaustive sweep of Rule 1' through
      W = 40 (FM-SEC152: 81.8M cases, zero failures) and adversarial
      search to W = 140.
    - Update 05:00 (FM-MECH143, FM-MECH144): lists with exactly two odd
      labels reduce to fewer factors (two-odd fusion), so the open part is
      pair-free lists with 0 or `>= 4` odd labels.  (D) and FM3 also hold on
      a uniform many-ones region.  (D) is exhaustive on the residual through
      W = 40 (FM-SEC164).
    - Update 05:30 (FM-SEC166): FLIP DESCENT.  If some double sign flip
      does not increase `Phi`, an exact positive identity reduces the list
      to lists with one fewer factor.  Through `W = 40` only 5,430 of
      21,364,508 residual cases lack one, and on all of them the fixed
      TopPair removal is monotone.  New single conjecture (FT): flip
      descent or TopPair.  It survives every known killer of the earlier
      rules.
    - No counterexample has been found anywhere.
  - **Checkpoint 2026-10-01 15:00 (main agent).**  FM3 full cone: NOT
    finished.
    - Proved and independently checked:
      - labels `<= 2` (FM-MECH41); labels `<= 3` (FM-MECH47, FM-CHK56);
        labels `<= 4` (FM-MECH48, FM-CHK57);
      - the many-factor theorems (FM-MECH49, FM-CHK59; constant weighted
        cutoff `T >= 2^21`, ADV-2, FM-CHK61);
      - (E) at `q = 2`, both signs, all labels (Q2+, FM-MECH50, FM-CHK58);
      - one arbitrary label with `b <= 2` copies of `hat S_2` (OL, Q2+, B2:
        FM-MECH51, FM-CHK60).
    - Proved, rerun by the main agent, no independent check yet:
      - (since checked by FM-CHK62: the `H = 1` layer at distance 3 (ADV-1)
        and distances 4, 5 (FM-MECH53); `min(a,e) <= 1` for every `b`, and
        the quadratic cutoff (FM-MECH52));
      - the (E) regions of FM-MECH55, 58, 59, 60 (balanced central
        anchors; far anchors; C1, C2) with outer propagation (FM-MECH57);
        FM-MECH55/56/57 are now independently checked (FM-CHK63), and
        FM-MECH58/59/60 too (FM-CHK64; the FM-MECH59 terminal step was
        repaired by an explicit identity).
    - Open:
      - (a) the `H = 1` layer at distance `>= 6` with `b >= 3`,
        `a, e >= 2` (per-row Gram certificates parked, FM-SEC136/137);
      - (b) (E) for `q >= 3` on the three ranges R1..R3 of FM-MECH60
        (running: FM-MECH61);
      - (c) two labels on a background containing `hat S_2` (`b >= 1`);
      - (d) few large labels on a background of many small labels `>= 3`;
      - (e) the labels `<= 5` box (227.3B entries; a fast port is running
        as FM-SEC139, to be run by the main agent with checkpoints).
    - The exact coverage atlas (FM-SEC138) puts the smallest uncovered
      lists at label sum 18, all containing a 5 together with a 3 or 4.
    - No counterexample has been found anywhere.
  - Finite-range rows.  The table rows restricted to finite `(r, a)` or
    finite `(a, e)` cover all labels there, since larger labels fall
    outside the support.  They are equivalent in strength to exact
    evaluation.  Their added value is that every case lies in a uniform
    criterion, which is evidence that the union of criteria is
    exhaustive.  The uniform results of item (39) are those stated for
    all parameters: EQ1, OL, LD, G0B, G0E, G0E4, the centered sector,
    and the metric and energy-drop regions.
  - Open: unbounded numbers of factors in every sector.  For T1, the
    product-class certificate fails at seven factors (bounded degree);
    slicing extends it slice by slice.
  - Q3 at level `r` needs exactly `2r` symmetric powers from the minus
    labels, plus `hat S` factors and `h_1`'s.  So the table covers:
    - Q3 words at `r = 1` in three cases:
      - one plus label `>= 2` and any number of plus 1s (FM53);
      - two plus labels `>= 2` and at most one plus 1 (S2-2, S2-3);
      - three plus labels `>= 2` and no plus 1 (S3-2);
    - Q3 words at `r = 2` with no plus label `>= 2` and at most two
      plus 1s, or one plus label `>= 2` and no plus 1;
    - no Q3 words at `r = 3` (they have `>= 6` factors).
- **Conjectured, with strong evidence.**
  - Net-cell GSI (item (38)): slot maps with the net cells of
    `prod hat S (x) Q_r`, with composite edges.  Distance 1 suffices for
    labels `(2,2)`, `(p,2)` and `(2,2,2)`, `(3,2,2)` at `r = 1`.
    Distance 2 suffices for `(3,3)`, `(p,1)`, `(4,3)`, `(5,3)`.  `(4,4)`
    at `r = 1` needs length 3.  It covers the `hat S` sectors where the
    two-layer lift fails.
  - GSI (item (38)): generic slot injection, H-only.  Checked:
    - `r = 2`: all words to size 12 except `(1^12)`, plus the size-14
      scan;
    - `r = 3`: size `<= 8`;
    - `r = 4`: 35 words.
    It implies FM3 in the H-only sector.
  - qFM3_2 (graded `r = 2` H-only): all 910 partitions with
    `|kappa| <= 18`.  Few-part screens (main agent), all with no negative
    coefficient:
    - at most 6 parts, `|kappa| <= 26`: 2,470 partitions;
    - at most 5 parts, `|kappa| <= 34`: 4,690 partitions;
    - at most 4 parts, `|kappa| <= 58`: 14,534 partitions.
  - (WITHDRAWN, see item (17)) Even invariant cohomology of the generic
    `r = 2` slot complex.  It holds for `|kappa| <= 8` but fails at
    `(4,1^6)`.
  - H-cone positivity at every integer `r`: all 1,597 partitions with
    `|kappa| <= 18` at every `1 <= r <= 10`, with no negative value
    (`verify_phir_hcone.py 18 1,...,10`).
  - Conjecture T (item (37)), the transport form of FM3: exact max-flow
    feasible on 86,376 H-only pairs and 19,788 + 13,120 `hat S` pairs.
  - Full FM3 cone (H-only and `hat S` words up to five plus factors),
    total degree `<= 22`, `r <= 20`: 285,020 (word, level) pairs, no
    negative value (item (35), closed-form screen validated against exact
    moments).


- **Killed.**
  - Analytic interpolation in `r`: non-integer levels fail.
  - A binomial expansion in `r`.
  - Factor-wise Koszul actions at `r = 3`: every orientation pattern at
    `(1,1,1,1)` has odd cohomology.
  - Graded fusion or split refinements at `r = 3`.
  - Graded level steps (linear programming).
  - A graded one-`hat S` refinement at `r = 2`: charge Kostka--Foulkes
    times ungraded `hat S_p` has negative coefficients for `p = 1..5`.
  - Local charge-shift injections for qFM3_2.
  - Row-exactness of `X` alone (item (6)).
  - Purity of the fusion Koszul complex (item (10)).
  - Square-zero and conjugated coupled actions at `r = 3` (item (11)).
  - Direct identification of `Phi` with the standard known-positive
    q-analogues (item (15)).
  - Evenness of any all-`+` slot complex at `r = 2` (item (17)), and the
    `A`-degree refinement at `r = 2` (size 10) and `r = 3`.
  - Cyclage-path matchings for qFM3_2 (item (19)).
  - Coefficientwise and log-concavity cones in Proposition 13 coordinates
    (item (24)).
  - Positive `h`/`e`/`p` expansions of `F_r` (item (25)).
  - Per-central-character positivity in the FM16 `Gr_2(C^4)` route
    (item (28)).
  - Local block matchings for the bigraded Euler characteristic (item
    (29)), and safe irreducible functionals for a Pieri induction (item
    (30)).
  - Integration-by-parts (Stein) induction on the level (item (32)).
  - Tail positivity of the central-character refinement as a uniform
    mechanism (item (33)): true for `r <= 5` in the tested range, false
    from `r = 6`.
  - Conjecture LP (item (34)): false at the hook `(14,1^14)`, which is
    negative on `(1.605, 1.784)`.  So no real-level positivity route
    exists.
  - Heckman--Opdam expansion positivity at `k = (1, r)` (item (35)):
    `h_2` already has a negative coefficient at `r = 2`.
  - Item (38), 2026-09-29:
    - per-`(m,n)` positivity of the `SU(2) x SU(2)` array of
      `(x-y)^(2r) h_kappa`;
    - per-shape payment in the level-1 two-label sector (FM-SEC2), and
      the per-shape `h`-expansion condition;
    - direct merge/slide multigraph injections for T1, at `(1^6)`
      (FM-PF1);
    - net-cell GSI with edges of bounded length for all labels at
      `r = 1`.
      - Witness: labels `(q+1, q+1)`, `kappa = (q, q)`.  The net cell
        `-2 V_(q,q)` has `m = 1`.  Within `q - 1` steps its only positive
        cell is `V_(q+1,q+1)`, and `m = 0` there; every other positive
        cell is `>= q` steps away.
      - Checked at `q = 3`: length 2 gives rank 3/5, length 3 gives
        5/5.

**(39) Absorbing the `h_1` factors of the consumer (FM-MECH1,
astra_medium_ceres; main-agent check `mech/n0check.py`, `mech/n0scan.py`).**
- Scripts cited as `mech/...`, `inj/...` and `xmn/...` in this item are
  in `character_ring_iter/fm39/` (see its README).
- **Consumer knob.**  Corollary FM2 needs, at level `r`, only words with at
  most `2r` general factors `h_k` (`k >= 2`), plus any number of
  `h_1 = S_1` and `hat S_p` factors.  FM3's unbounded number of general
  `h_k` is not used by the consumer.
- **Reduction (confirmed).**  Write a core word as
  `w = s^eps sum c_ab s^(2a) d^(2b)`, with `s = x+y`, `d = x-y`.  Then
  `phi_r(h_1^(2m+eps) w) = sum c_ab M(r+b, m+a+eps)`, where
  `M(r,m) = phi_r(h_1^(2m))` is the all-ones value.  So
  `phi_r(h_1^(2m+eps) w)/M(r,m)` is rational in `m`.  Checked from
  Catalan moments, without the closed form of `M`.
- **KILL (N0, zero-shift Newton positivity).**  With `P = D(m) x (that
  ratio)` and `D = (m+r+2)^(L up) (m+r+3)^(L up)`, all Newton coefficients of
  `P` are `>= 0` for Astra's four cores.  But not in general.
  - Screen: up to `2r` factors `h_k` (`2 <= k <= 7`) times up to three
    `hat S_p` (`2 <= p <= 5`).  Degree `<= 12` at `r = 1, 2`, `<= 10` at
    `r = 3, 4`.
  - Failures: 0/221, 1/296, 5/149, 13/149 at `r = 1..4`.
  - First failure: `hat S_4` at `r = 2`, with coefficients
    `[120, 1560, 1464, -120, 240]`.  Others: `hat S_2^2` and `h_2 hat S_2^3`
    at `r = 3`; `hat S_4^2` and `h_2 hat S_2^2` at `r = 4`.
  - `phi_r(h_1^n w) >= 0` for `n < 80` on every failing core.  Only the
    certificate fails, not positivity.
  - The knob violated is uniformity over cores at shift 0.  A
    core-dependent shift survives, but it certifies one core at a time,
    and the cores are still an infinite family.
- The other two FM-MECH1 suggestions refine routes already recorded: cut
  compression (the band-minimizer conjecture, whose interval positivity
  still contains T1) and a surviving determinant coefficient (the generic
  rank in GSI).  FM-MECH2 (the same thread at max effort) is assigned to
  label-uniform mechanisms.
- **Theorem R3-S1 (FM-MECH2, astra_max_ceres; verified by the main agent,
  `mech/kraw_check.py`, `kraw_sym.py`, `kraw_bdry.py`).**
  `phi_3(hat S_p h_1^a) >= 0` for every `p >= 2` and `a >= 0`.  This
  closes the paused one-label residual (FM27, FM33, FM36, FM37) at
  `t = 2r = 6`.
  - Put `N = a+6`, `k = (N+p)/2`, `n = (N-p)/2`.  The value is 0 off
    parity or for `p > N`.  By FM27 it is `D_k - D_(k+1)`, with
    `c_j = [z^j](1+z)^a (1-z)^6`.
  - Identity: `D_k - D_(k+1) = C(N,k)^2 (p+1)(N-5)(N-4) R_6(N, p(p+2)) /
    [(n+1)(k+1)^2(k+2)(N(N-1)...(N-5))^2]`.
    - Checked as an exact rational identity in `(a,k)` in the interior.
    - Checked as a polynomial identity in `N` for each boundary family
      `n = 0..7`.
    - It matches direct Catalan-moment values at 888 points
      (`8 <= N <= 60`).
  - Positivity: for `N >= 8`, with `z = 8/N` and `u = p(p+2)/N`,
    `R_6(N, Nu)/N^6 = sum_(i=0..4) b_i(u) C(4,i) z^i (1-z)^(4-i)`.  Each
    `b_i` is a degree-6 polynomial with positive constant term and no real
    root in `[0, oo)`.  The integer coefficient vectors are:
    - `b_0`: `1575, -3150, 4095, -1740, 345, -30, 1`;
    - `16 b_1`: `20475, -44730, 58710, -25800, 5235, -470, 16`;
    - `96 b_2`: `96075, -232830, 312525, -142620, 29752, -2760, 96`;
    - `256 b_3`: `193725, -527670, 732030, -348432, 75056, -7200, 256`;
    - `128 b_4`: `70875, -220050, 318420, -158784, 35456, -3520, 128`.
  - Base cases: `N = 6` gives `84, 20, 1` (`p = 2, 4, 6`); `N = 7` gives
    `40, 15, 1` (`p = 3, 5, 7`).
  - Scope: one fixed level.  The Krawtchouk order and the number of
    coefficient polynomials grow with `r`, so this is not a uniform-in-`r`
    mechanism.
- **Theorem R3-TWO (FM-MECH2, astra_max_ceres; verified by the main agent,
  `mech/r3tail_check.py`, `r3tail_check2.py`).**
  `phi_3(h_s h_t h_1^a) >= 0` for all `s, t, a >= 0`.  This closes the
  last region `n >= 161`, `0 <= d < n` of the `r = 3` two-power family
  (FM-T1-TR3-CERT), whose other regions are proved above.
  - With `A = SQ` and `C_1, C_2, C_3` as in FM-T1-TR3-CERT, put
    `J = 243 A C_2^2 - 270 A C_2 C_3 - 25 A C_3^2 + 75 A^2 C_2 + 80 C_3^3`
    and `H(x) = A x^3 + 5A x^2 + 20 C_3 x + 60 C_2`.
  - Degree-5 Taylor: `e^x - P(x) >= x C_1/(2A) + x^2 H(x)/(120A)` for
    `x >= 0` (exact identity against `T_5`), and
    `disc(H) = -400 A J` (exact).
  - So `P <= e^x` on `x >= 0` whenever `C_1, C_2 >= 0` and
    (`C_3 >= 0` or `J >= 0`).
  - These hold for `n >= 161`, `d < n`:
    - `n >= 169`, `d >= 4 sqrt(n)`: under `n = (13+v)^2`,
      `d = 4(13+v)+e`, every coefficient of `C_1, C_2, C_3` is positive
      (minimum 1).
    - `n >= 169`, `d <= 4 sqrt(n)`: with `z = d/sqrt(n)` and
      `h = 1/sqrt(n)`, the Bernstein coefficients on `h` in `[0, 1/13]`
      are all `>= 0`, with no subdivision:
      - `z` in `[0,2]`: `h^9 C_1`, `h^10 C_2`, `h^30 J`;
      - `z` in `[2,4]`: `h^9 C_1`, `h^10 C_2`, `h^10 C_3`.
    - `161 <= n <= 168`: 1,316 exact checks.
  - This is the `r = 3` two-power slice that was stopped at the user's
    direction (direction note above).  It is now complete, but it is still
    a per-family result.
- **One-label family at every level: a uniform reduction (FM-MECH3,
  astra_max_ceres; verified by the main agent, `mech/gram_repro.py`,
  `gram_link.py`, `gram_ext.py`).**
  - *Factorization, every `t`.*  With `p = 2k-N`, `n = N-k`,
    `T = p(p+2)`, and Krawtchouk polynomials `K_0 = 1`, `K_1 = p`,
    `K_(j+1) = p K_j - j(N-j+1) K_(j-1)`:

        D_k - D_(k+1) = C(N,k)^2 (p+1)(N-t+1)(N-t+2) R_t(N,T)
                        / [(n+1)(k+1)^2(k+2) ((N)_t)^2],

    where `R_t(N, p(p+2)) = A_t/[(p+1)(N-t+1)(N-t+2)]` and
    `A_t = (n+1)(k+1)^2(p+2) f^2 + p n^2 (k+2) g^2 - (N-2t)(p+1) n (k+1) f g`,
    with `f = K_t(p;N)` and `g = K_t(p+2;N)`.  `R_t` is monic of degree `t`
    in `T`.
    - Checked: 513 coefficient comparisons (agent code, rerun).
    - Checked against direct Catalan-moment values of `phi_r` for both
      parities (`hat S_q` at `t = 2r`, `Sym^s` at `t = 2r-1`), `t = 4..10`:
      616 values, no mismatch.  This also confirms FM27 up to `r = 5`.
  - *Folding.*  Swapping `a` and `t` leaves `D_k` unchanged, so it suffices
    to treat `N >= 2t` (with `t` replaced by `min(a,t)`).  `min <= 3` is
    FM33.
  - *Prescribed sum of squares.*  With `P_j(X) = K_(t-2j)(X; N+2)`, the
    triangular expansion
    `R_t(N,X^2) = sum alpha_j P_j^2 + 2 sum beta_j P_j P_(j+1)` is unique.
    - With `rho_j = (t-2j)(N-t+2)` and
      `delta_j = alpha_j + rho_(j-1) beta_(j-1) + beta_j/rho_j`, it becomes
      `R_t = sum delta_j P_j^2 + sum (-beta_j/rho_j)(P_j - rho_j P_(j+1))^2`.
      This is an exact identity.
    - So `R_t >= 0` for all `T >= 0` if `beta_j <= 0` and `delta_j >= 0`.
  - *Proved for all `t`:* rows `j = 0, 1`.
    - `alpha_(t,1) = t(t-1)[(N-t+2)(N-t+3) - 2(t+1)/3]`.
    - `beta_(t,1) = -2t(t-1)(t-2)(t-3)(5N-2t+18)/15`.
    - `E_(t,1)(2t+v)/(t(t-1)(t-2))` has all coefficients positive in
      `(t, v)`.
  - *Checked:* every row, all `t <= 32` (256 sign polynomials, all
    coefficients `>= 0` after `N = 2t+v`).
  - *Rows 2, 3, 4 proved for all `t` (FM-MECH4, astra_max_ceres; rerun by
    the main agent, `mech/gramrows.py`).*
    - Normalize `alpha_(t,j) = (t)_(2j) A_j(t,h)` and
      `beta_(t,j) = (t)_(2j+2) B_j(t,h)`, with `h = N-t+2`.
    - Put `F_j = h A_j + (t-2j+2) h^2 B_(j-1) + (t-2j-1) B_j`.  The sign
      rows are then
      - `E_(t,j) = (t)_(2j)(t-2j) F_j` for interior rows;
      - `E_(t,j) = (t)_(2j) F_j/h` at `t = 2j+1`;
      - `E_(t,j) = (t)_(2j) Z_j(N-2t)` at `t = 2j`.
    - Three polynomials carry the signs:
      - `O_j(u,v) = -B_j(2j+2+u, 2j+4+u+v)`;
      - `P_j(u,v) = F_j(2j+1+u, 2j+3+u+v)`;
      - `Z_j(v) = A_j(2j, 2j+2+v) + 2(2j+2+v) B_(j-1)(2j, 2j+2+v)`.
    - For `j = 1..7` every coefficient of `O_j`, `P_j`, `Z_j` is strictly
      positive, symbolically in `t`.  Rows 1--4 are FM-MECH4; rows 5--7 are
      the main agent's run of the same constructor (`mech/gramrows7.log`).
      - The coefficient counts are `j(2j+1)`, `(j+1)(2j+3)` and `2j+1`.
      - The minima are `2j/3`, 1 and 1.
    - So the sign lemma holds for rows `j <= 7` at every `t`, and for every
      row at `t <= 32`.  Open: rows `j >= 8` with `t >= 33`.
    - `Z_j(v)` is monic of degree `2j`, and its next coefficient is
      `3j(2j+1)` for `j <= 7`.
    - FM-MECH5 (astra_max_ceres; reproducer rerun by the main agent):
      - A closed binomial inverse for the Gram coefficients in the basis
        `K_d(X; N+2)^2`, with a two-term row recurrence that replaces the
        chain sums (98 + 97 exact checks).
      - An explicit product-weighted formula for `Z_j`, with a positive
        `3F2` kernel.  The kernel alone does not give positivity, since the
        functional is not positive on squares.
      - A proof that `v + 4j + 2` divides `Z_j(v)` for every `j`.
      - Signs still open.
    - Closed finite formulas exist for every row.  `r_l = [T^(t-l)] R_t` is
      a coefficient extraction from an explicit series, and
      `alpha`, `beta` are chain sums (the inverse of a unit triangular
      matrix).  These sums alternate, so they do not give the signs
      directly.
    - KILL (one transfer): `R_t - K_t^2 - t(t-1)(N-t+1)(N-t+2) R_(t-2)` is
      not Gram-positive in the inherited basis.  At `N = 2t` its row-2
      diagonal is `-2(114t^3 - 3929t^2 - 10233t - 3250)/315`, which is
      negative at `t = 64`.  The knob violated is positivity of that
      remainder matrix, which the consumer does not need.  The sign lemma
      itself is not affected.
  - *Hermite limit:*
    `N^(-t) R_t(N,Nu) -> sum_j t!/(t-2j)! He_(t-2j)(sqrt u)^2`.
  - **Open (first gap).**  The Gram-row sign lemma for rows `j >= 2`,
    uniformly in `t`: `-beta_(t,j)(2t+v) >= 0` and `E_(t,j)(2t+v) >= 0`.
    Proving it would close the one-label sectors (`hat S_q h_1^a` and
    `Sym^s h_1^a`) at every level `r`.
- **Theorem OL (the one-label sector at every level; FM-MECH6,
  astra_max_ceres; verified by the main agent).**  For all integers
  `a, t >= 0`, `N = a+t`, `c_j = [z^j](1+z)^a(1-z)^t` and
  `D_k = c_k^2 - c_(k-1) c_(k+1)`, one has `D_k >= D_(k+1)` whenever
  `2k > N`.  With FM27, whose Catalan identity was accepted with a proof in
  FM-CHK9:

      phi_r(hat S_q h_1^a) >= 0   and   phi_r(h_s h_1^a) >= 0

  for every level `r` and all `q, s, a` (`t = 2r` and `t = 2r-1`).  So Q3
  holds for every word whose labels are all 1 except one arbitrary label,
  at every level.  This closes the residual branches `R_-`, `R_+`, which
  had been paused after three attempts.
  - *Proof outline.*  Fold to `N >= 2t` (swapping `a` and `t` leaves `D_k`
    unchanged).  `t <= 3` is FM33, and `p = 2k-N >= (N-2)/3` is FM27b.  On
    the rest (`1 <= p < (N-2)/3`) put `s = p+1`, `n = N-k`,
    `d = N-2t`.
    1. `(n+1)(k+2)(D_k - D_(k+1)) = c_k^2 Q(c_(k+1)/c_k)`, with
       `Q(r) = p(k+2) r^2 - s d r + (p+2)(n+1)`.  If
       `F = s^2 d^2 - 4p(p+2)(n+1)(k+2) <= 0`, `Q >= 0`: the FM33 band.
    2. *Outer region*, `d^2 >= 4(n-1)(k+2)`.
       - `c_(k+1)/c_k = n K_(n-1)(d;N)/K_n(d;N)`.
       - The continued fraction `u_j = 1/(d - a_(j-1) u_(j-1))`, with
         `a_j = j(N-j+1) <= (n-1)(k+2)`, gives
         `0 < r <= R = 2n/(d + sqrt(d^2 - 4(n-1)(k+2)))`, and `R` lies left
         of the vertex of `Q`.
       - `(n-1) Q(R) >= (2n-1)p + 2(n-1) > 0`.
    3. *Reciprocal roots.*  The zeros of `K_t(x;N)` are real, simple and
       symmetric (Jacobi matrix).  For `N >= 2t`,
       `sum 1/lambda_i^2 <= t/(2(N-t+2))` for even `t`, and
       `<= (t-1)/(4(N-t+3))` for odd `t` (nonzero roots).  The proof uses
       the three-term recurrence at `x = 0` and a weighted-mean (Chebyshev
       sum) bound.
    4. *Central region*, `N >= 3ts^2`.  Then `p, p+2` lie below the first
       positive zero, and by duality
       `c_(k+1)/c_k = (n/(k+1)) K_t(p+2;N)/K_t(p;N)`.
       - Even `t`: `0 < r <= 1 <=` vertex, so `Q(r) >= Q(1) = 2(t+1)s`.
       - Odd `t`: `r >= U(1 - bs/N)` with `U = (s+1)/(s-1)` and
         `b = (4t+2)/3`, which lies right of the vertex.  The tangent at
         `U` gives `Q(r) >= Us (4t^2 + 26t - 15)/(18t) > 0`.
    5. *Coverage.*  If `F > 0` and `d^2 < 4(n-1)(k+2)`, then `N >= 4s^2`.
       `F` is an upward quadratic in `N` that is negative at
       `max(2t, 4s^2)` and at `3ts^2`, so `N > 3ts^2`.
  - *Main-agent verification* (`mech/onelabel_repro.py`,
    `mech/onelabel_indep.py`):
    - Every step was re-derived by hand, including the recurrences at
      `x = 0`, the Chebyshev bound, the vertex comparisons and the
      coverage identity `F = s^2(d^2 - 4(n-1)(k+2)) + 4(k+2)(n+1-2s^2)`.
    - The agent's vetting code reruns exactly.
    - Independent exact checks:
      - coverage at 42,614,378 middle-range points (`t <= 60`,
        `N <= 3000`), none uncovered;
      - the reciprocal-root bounds at 30,877 pairs (`t <= 80`,
        `2t <= N <= 2t+400`), no failure;
      - the theorem together with each lemma's intermediate claim at
        350,397 points (`t <= 30`, `N <= 400`; band 343,444, outer 4,741,
        central 2,212), no failure.
    - The factorization and direct `phi_r` values agree for
      `t = 4..10` (616 values).
  - The Gram-row sign lemma above is not needed.  It stays open as a
    stronger statement.
- **Two labels at every level: one inequality (FM-MECH7,
  astra_max_ceres; verified by the main agent, `mech/twolabel_repro.py`,
  `mech/twolabel_indep.py`).**
  - *Reduction (exact).*  Let `c_k = [z^k](1+z)^a(1-z)^e`, `N = a+e`,
    `B_k = c_(k-1) + c_(k+1)`, `D_k = c_k^2 - c_(k-1) c_(k+1)`.  For
    labels `p >= q` put `j = (N+p-q)/2` and `i = (N+p+q)/2 + 1`.  Then
    - the cross kernel is `E[(x-y)^e (x+y)^a U_p(x) U_q(y)] = W_ij`, with
      `W_ij = B_i c_j - c_i B_j`;
    - its row interval is `T_ij = D_j - D_i`.
  - The two-label words then evaluate as follows (labels sorted):
    - `h_s h_t h_1^a` at `e = 2r-2` (labels `s+1, t+1`) is `T - W`;
    - `hat S_p hat S_q h_1^a` at `e = 2r` is `T + W`;
    - `h_s hat S_q h_1^a` at `e = 2r-1` (labels `s+1, q`) is `T +- W`.
  - Each of these is a 2x2 determinant of shifted `c`'s.  So the whole
    two-label sector at every level is the single inequality

        (E)   D_j - D_i >= |B_i c_j - c_i B_j|,   N/2 <= j < i.

    - It is invariant under the fold `a <-> e`.
    - In Krawtchouk form it reads
      `G_j - rho^2 G_i >= rho |H_i f_j - f_i H_j|`, with
      `rho = C(N,i)/C(N,j)`.
  - *Proved at every level:* `|a - e| <= 1`, by exact factorizations.
    - At `a = e`, the slack is `(A-B)(A+C)` or `(B-C)(A+C)` in central
      binomials.
    - At `a = e+1`, it is `(x-y)(x+y+A+C)` or `(x+y)(x-y+A-C)`.
    - Also `min(a,e) <= 2`, by folding to FM42, FM49 and FM54.
    - So these word families hold at every `r` and for all labels:
      - `h_s h_t h_1^a` for `a <= 2` or `a` in `{2r-3, 2r-2, 2r-1}`;
      - `hat S_p hat S_q h_1^a` for `a <= 2` or `a` in `{2r-1, 2r, 2r+1}`;
      - `h_s hat S_q h_1^a` for `a <= 2` or `a` in `{2r-2, 2r-1, 2r}`.
  - *Checks.*
    - The reduction agrees with the main agent's `phi` in 2,700 cases
      (three patterns, `r <= 4`).
    - The strip factorizations were re-derived by hand.
    - `(E)` holds at 38,826 pairs on the strip (`e <= 40`) and at 71,527
      pairs for all `a, e <= 24`, including the open region.
  - KILL (monotone extension): `T_ij - W_ij` is not monotone in `i`.  At
    `a = 8`, `e = 6`: `phi_4(h_6 h_2 h_1^8) = 1064 > 888 = phi_4(h_7 h_3 h_1^8)`.
    The knob violated is the monotonicity, which the consumer does not
    need.
  - **Open:** `(E)` for `|a - e| >= 2`, `a, e >= 3`.
  - **Root-free intervals (FM-MECH8, astra_max_ceres; verified by the main
    agent, `mech/rootfree_repro.py`, `mech/rootfree_indep.py`).**
    - Fold to `a >= e >= 3`.  Put `n = N+2`, `X_k = 2k-N`, `z_k = X_k^2`,
      `b_k = C(n,k+1)`, `nu_l = l! (n)_l`, and
      `eta = e!/[4(n)_(e+2)]`.
      Here `(n)_l = n(n-1)...(n-l+1)` is the falling factorial (FM-CHK26:
      a rising reading breaks the identity, e.g. at `e = 1`, `a = 3`,
      `(j,i) = (3,4)`).
    - Define the positive kernel
      `Gamma(u,v) = b_u b_v sum_(l <= e, l = e mod 2) K_l(X_u;n) K_l(X_v;n)/nu_l`.
    - Christoffel--Darboux in parity form gives
      - `W_ij = eta (z_i - z_j) Gamma(i,j)`;
      - `D_j - D_i = eta sum_(k=j)^(i-1) (z_(k+1) - z_k) Gamma(k,k+1)`.
    - Both identities were checked with the main agent's own code on 8,200
      pairs.
    - Weighted chord lemma: if `log g_k` is concave in `z_k`, then
      `sum (z_(k+1)-z_k) g_k g_(k+1) >= (z_i - z_j) g_j g_i`.  This is
      weighted AM--GM, with `sum w_k (theta_k + theta_(k+1)) = 1` exactly.
    - `log b_k` is concave in `z` (exact slope
      `-atanh((X+1)/(n+1))/(2(X+1))`).  `log|K_l(sqrt z)|` is concave
      between roots.
    - **Theorem RF.**  If `(X_j, X_i)` contains no root of any `K_l(X;n)`
      with `l <= e`, `l = e mod 2`, then `D_j - D_i >= |W_ij|`, for both
      signs.
    - Explicit regions follow, using Theorem OL's reciprocal-root bounds
      and a Gershgorin bound on the spectral radius:
      - central: `Y^2 <= 2(n-m+2)/m` (`m` even) or
        `Y^2 <= 4(n-m+3)/(m-1)` (`m` odd);
      - outer: `X^2 >= 4(m-1)(n-m+2)`.
    - **Corollary.**  For every `r >= 2`,
      `a >= (r-1)(s+t+4)^2 - 4` implies `phi_r(h_s h_t h_1^a) >= 0`.  The
      same regions apply to the plus-plus and mixed patterns.
    - Termwise fails across roots.  At `(a,e,j,i) = (8,6,8,13)` the
      degree-4 term is `-432/11`, while the sum is `T - W = 1176`.  So only
      the sum over `l` can be positive there.
    - Open: intervals `(X_j, X_i)` that contain a root.
  - **Positive-sum form (FM-MECH9, astra_max_ceres; reproducer rerun and
    scan by the main agent, `mech/pe_repro.py`, `mech/pe_scan.py`).**
    - Setup.  On the squared lattice `z_k = X_k^2` (right half, including
      `k = N+1`), put `mu_k = (2 - 1_(X_k=0)) b_k X_k^(2 eps)`.
    - By Cauchy--Binet, the kernel is a positive Vandermonde ensemble:
      `Gamma(u,v) = (2^n/Z) g_u g_v H_q((z-z_u)(z-z_v))`, with
      `H_q(f) = sum_(|S|=q) mu(S) Delta(S)^2 prod_S f`.
    - Adjacent nodes give only nonnegative summands.  This is a short
      lattice proof of Theorem OL's adjacent case, across roots.
    - Stronger target (PE): replace `|H_q(...)|` for the chord by
      `H_q(|...|)`.  (PE) implies `(E)` for both signs, and it is
      invariant under the fold (particle--hole duality of the binomial
      lattice).  At even `e` it reads
      `sum_(k=j)^(i-1) Lambda_(k,k+1) >= Lambda_(j,i)`, where
      `Lambda_uv` is a sum of positive minors.
    - Not provable subset by subset.  One subset has every adjacent product
      0 and a positive chord.
    - Not true for arbitrary log-concave grids.  On `z = (0,1,10,19,20)`
      with unit weights: 760 against 2000.
    - Exact scan: (PE) holds at all 22,944 pairs with `e = 3..10`,
      `a >= e+2`, `N <= 40`.
    - Open: (PE), or `(E)`, across roots.
  - **Two-step splitting (FM-MECH10, astra_max_ceres; reproducer rerun,
    witness value `phi_4(h_3 h_2 h_1^17) = 5913807` recomputed by the main
    agent).**
    - (PE) holds for every interval of two lattice steps, at every level
      and for both parities.
      - The proof groups the subsets by which of the three nodes they
        contain.
      - Each group is then a convexity statement for `g`, or for the dual
        weight `g/(mu |F'|)`, and both are log-concave on the binomial
        squared lattice.
    - Exact deletion recursion: `P = P^- + mu_h P^+ + delta_h`, with
      `delta_h >= 0`.
    - KILL (branchwise positivity), H-only: at `(a,e,j,h,i) =
      (17,6,12,14,16)`, i.e. `h_3 h_2 h_1^17` at `r = 4`, `mu_h P^+ < 0`
      while `P > 0`.  The knob violated is separate positivity of each
      branch, which the consumer does not need.
    - KILL (local-to-global): on a flat log-concave grid all triple slacks
      are positive, but the path is 128 against a chord of 192.  So the
      proof must use the lattice.
    - **Circuit breaker (2026-09-30).**  FM-MECH9 and FM-MECH10 were two
      consecutive auxiliary batches on the prerequisite (PE)/(E), and
      neither discharged a Q3 obligation.  Further auxiliary assignment on
      this detour needs explicit user authorization.
    - Last movement: Theorem OL, the strip and Theorem RF each prove Q3
      directly for infinite word families.
  - **KILL (norm-chord route; main agent, `mech/normchord.py`).**  The
    Cauchy--Schwarz strengthening
    `T_ij >= eta (z_i - z_j) sqrt(Gamma(i,i) Gamma(j,j))` fails on 5,465 of
    65,455 pairs (`e = 3..12`); worst ratio 0.23 at `(a,e,j,i) =
    (15,12,14,16)`.  The knob violated is replacing the kernel by its
    norms, which the consumer does not need.  (PE) sits strictly between.
  - **Falsification screens after the circuit breaker (main agent, exact).**
    - `(E)` itself: 15,258,258 pairs with `e <= 60`, `a <= 150`; no
      failure (`mech/escreen.py`).  The run to `e <= 120`, `a <= 400` is
      recorded below.
    - FM3 itself, through the kernel reduction
      `phi_r(w h_1^a) = (1/2) sum_(p,q) F_pq(w) W(p,q; a, 2r)`.
      - The reduction was validated against direct `phi` in 224 cases
        (`mech/fm3kern.py`).
      - Over all cores of up to 5 labels from `h_2..h_6`,
        `hat S_2..hat S_6` (`r <= 20`, `a <= 100`), and up to 3 labels from
        `h_2..h_8`, `hat S_2..hat S_8` (`r <= 40`, `a <= 150`):
        5,084,040 exact evaluations, no negative value
        (`mech/fm3screen.py`).
  - **H-only consumer as a split sum.**  For `m <= 2r` general parts
    `A_i = kappa_i + 1`, absorb every `(x-y)`:
    `phi_r(h_kappa h_1^a) = (1/2) sum_(S subset [m]) (-1)^(m-|S|)
    <U_S (x) U_(S^c), W_(2r-m, a)>`, where `U_S = prod_(i in S) U_(A_i)`.
    - `m = 1` is Theorem OL.
    - `m = 2` is `(E)`.
    - The consumer's cap of `2r` general factors is exactly what makes
      `2r - m >= 0`.
  - **Theorem THREE (FM-MECH11, astra_max_ceres, direct construction;
    verified by the main agent, `mech/parity_check.py`,
    `mech/three_repro.py`).**  `phi_r(h_u h_v h_w) >= 0` for every level `r`
    and all `u, v, w >= 1`.  This closes the FM-THREE aggregate.
    - *Whole-word parity transformation.*
      - Let `t` be the number of odd parts and `e = 2r - m`.
      - Under `y -> -y` (duality (8)), `x-y` and `x+y` swap.  An odd part
        stays a difference `U(x) - U(y)`, and an even part `kappa` becomes
        `hat S_(kappa+1)`.
      - So for `a + t` even and `R = (a+t)/2 >= 1`:
        `phi_r(prod h_(kappa_i) h_1^a) = phi_R(prod_(odd) h_(kappa_i)
        prod_(even) hat S_(kappa_i + 1) h_1^e)`.
      - If `a + t` is odd, the value is 0.
      - If `a = t = 0`, the value is `(1/2) E[(x+y)^e prod hat S]`, a genuine
        character, so it is `>= 0`.  So every H-only consumer word with all
        parts even and no suffix is `>= 0`, for every `m <= 2r`.
    - *Level 1.*  `phi_1(h_u h_v h_w) = <h_u h_v h_w, 1>_Sp4 >= 0` directly, since
      the transformed suffix `2r-3` would be negative there (FM-CHK26
      repair).
    - *Three factors, `r >= 2`.*  Two odd parts and one even part map to level 1 with
      one `hat S`: `phi_r(h_u h_v h_c) = phi_1(h_u h_v hat S_(c+1)
      h_1^(2r-3))`, which is `>= 0` by FM53.  Every other parity pattern is
      0 or all-even.
    - *Checks.*
      - The transformation was checked by the main agent in 1,580 cases
        (`r <= 4`, up to 4 parts, `a <= 3`), with no mismatch.
      - `phi_r(h_u h_v h_w)` was recomputed for `r <= 8` and parts `<= 7`
        (672 values), none negative.
      - The agent verifier reruns, including an explicit positive
        tableau count for the two-odd case.
    - *First suffix:* `phi_r(h_s h_(p-1) h_(q-1) h_1)` becomes a two-`hat S`
      word at level 1.  Its shape-by-shape count has the negative family
      `C_pq(p, q-1, 1, 0) = -1`; for example
      `phi_4(h_2 h_3 h_4 h_1) = 91 - 20 = 71`.  So per-shape positivity
      fails there, not the consumer.
  - (E) screen, completed: every `e <= 120`, `a <= 400` (all index pairs,
    folded), no failure.
- **Theorem LS (large suffix, the whole FM3 cone; FM-MECH12,
  astra_max_ceres; verified by the main agent, `mech/radial_repro.py`,
  `mech/ls_check.py`).**
  - For every word `w`, a product of any `h_kappa` and `hat S_p`, every
    level `r` and every `a >= 0`, put
    `Lambda(w) = sum_(h) kappa(kappa+4)/5 + sum_(hat S) p(p+2)/3`.  Then
    `phi_r(w h_1^a) >= dim(w) Z_(r,a) [1 - (2r+3) Lambda(w)/(a+2r+4)]`,
    with `Z_(r,a) = (1/2) E[(x-y)^(2r) |x+y|^a]`, which equals
    `phi_r(h_1^a)` for even `a`.  This corrects `Z = phi_r(h_1^a)`, which is
    wrong for odd `a` (FM-MECH13).  The odd-`a` form was checked in 225
    cases.
  - The bound is stated for `deg(w) + a` even.  If `deg(w) + a` is odd,
    `phi_r(w h_1^a) = 0` by the sign reversal `(x,y) -> (-x,-y)`.  FM-CHK26
    rejected the version without this clause: `w = h_1`, `a = 0` gives
    `phi = 0` against `2/3`.  With the clause and `Z = Z^+ =
    int_(x+y>0) (x-y)^(2r) (x+y)^a`, the proof applies.
  - So `a + 2r + 4 >= (2r+3) Lambda(w)` implies `phi_r(w h_1^a) >= 0`.
    For each core and each level, only finitely many suffix lengths
    remain.
  - *Proof.*
    - `w` is a genuine `SU(2) x SU(2)` character.  Its torus weights have
      mean 0 and coordinate variance `V = Lambda/2`
      (`kappa(kappa+4)/10` per `h_kappa`, `p(p+2)/6` per `hat S_p`).
    - The bounds `1 - cos A cos B <= (1 - cos A) + (1 - cos B)` and
      `1 - cos(I theta) <= I^2 (1 - cos theta)` give
      `w/dim(w) >= 1 - V(4-s)/2`, with `s = x+y`.
    - On `s > 0` (sign reversal), the law of `s` under
      `(x-y)^(2r) s^a` has density `s^a (4-s)^(2r+2) g_r(s)`.  Here
      `g_r(s) = int_0^1 z^(2r) sqrt(1-z^2) sqrt((4+s)^2 - (4-s)^2 z^2) dz`
      is increasing.  Comparison with a Beta law then gives
      `E[4-s] <= 4(2r+3)/(a+2r+4)`.
  - *Checks.*
    - The density was re-derived by the main agent (substitute
      `d = (4-s)z`).
    - The quantitative bound holds exactly in 1,760 cases with `h` and
      `hat S` factors (`r <= 4`, `a < 40`).
    - The mean bound holds numerically, worst ratio 0.988.
- **Theorem LR (large `r + a`, H-only; FM-MECH12; verified).**  For
  `w = prod h_(kappa_i)` with `t` odd parts, `d = (sum kappa - t)/2`,
  `B = 4^(-t) prod C(kappa_i+3, 3)` and `K = 16 B sum_(j<=d) j^2`:
  `phi_r(w h_1^a) >= phi_r(h_1^(a+t)) [1 - K/(2r+a+t+6)]`.
  - *Proof.*
    - Write `w = (x+y)^t F`, using
      `h_(2j+1) = (x+y)[h_j + h_(j-1)](x^2-2, y^2-2)`.
    - `F >= 1` on the boundary of `[-2,2]^2`, by
      `h_(2j)(2,y) = sum (U_l + U_(l-1))^2` and
      `h_(2j+1)(2,y) = (2+y) sum U_l^2`.  Also `|F| <= B`.
    - Chebyshev coefficient bounds give `P_q(z) >= 1 - L(1-z)` in
      sector coordinates `x = 2 sqrt z`, `y = 2q sqrt z`.
    - Radial comparison with `Beta(N/2+1, 2)` gives
      `E[1-z | q] <= 4/(N+6)`.
  - The identities were checked symbolically for `j <= 4`, and the bound
    was rerun in 315 cases.
  - Consequence: for each H-only core, FM3 can fail at only finitely many
    `(r, a)`, namely `2r + a + t + 6 < K`.  But `K` grows like a product
    of dimensions, and the cores are unbounded.
- **FM-CHK26 (luna_max_venus, independent checker, own code) on item (39):**
  - OL: ACCEPT.  Every step was re-derived, with 61,900 recurrence
    checks, 84,647 quadratic-identity checks and 2,929 root-bound checks.
  - LR: ACCEPT (5,439 exact cases).
  - RF: repaired (falling-factorial convention made explicit; 2,747
    Christoffel--Darboux pairs and 195 root-free intervals checked).
  - THREE: repaired (`r = 1` done directly; 6,020 parity checks).
  - LS: the version without the parity clause was rejected.  It is
    repaired above, and the large-suffix positivity consequence holds.
  - LR4 was not in the checker's scope.
- **FM-CHK27 (luna_max_venus, own code):** LS (repaired) ACCEPT; LR4 ACCEPT; BW ACCEPT, with the clarifications above.
- **Theorem LR4 (additive quartic cutoff, every word; FM-MECH13,
  astra_max_ceres; verified by the main agent, `mech/lr4_repro.py`,
  `mech/lr4_check.py`).**
  - Let `w` be any product of `h_k` and `hat S_p`, `t` its number of
    odd labels, and `N = 2r + a + t`.  Put
    - `l_H(k) = 2k^2(k+1)(k+3)/3` for even `k`, `(k-1)^2(k+2)(k+3)/6` for
      odd `k`;
    - `l_S(p) = 2p^2` for even `p`, `l_H(p-1)` for odd `p`;
    - `L(w) = sum l`.
  - Then `N + 6 >= 7 L(w)` implies
    `phi_r(w h_1^a) >= phi_r(h_1^(a+t))/25 > 0` (with `a+t` even; for `a+t`
    odd the value is 0).  Also `7L <= (14/3)(D^4 + 4D^3 + 3D^2)`, `D` the
    label sum, independently of the number of factors.
  - *Proof.*
    - Each factor's radial polynomial `P_q` satisfies
      `|P_q(z)/P_q(1) - 1| <= l (1-z)`.
      - This combines Markov's inequality `||P'|| <= 2j^2 ||P||` with
        dimension bounds.
      - It also uses boundary lower bounds `P_q(1) >= (j+1)/4` (or
        `(2j+1)/2`), from the pairing identity
        `U_l^2 + U_(l-1)^2 - y U_l U_(l-1) = 1`.
      - Odd `hat S` factors reduce via
        `hat S_(2j+1) = (x+y) h_(2j)(x,-y)`.
    - For the product,
      `G >= 2 - exp(L(1-z))`.  Beta moments give
      `E exp(L(1-z)) <= (1 - 2L/(N+6))^(-2) <= 49/25`.
    - Normalization (FM-CHK27).  Keep the boundary product
      `B(q) = prod P_i(1,q) >= 1` inside the expectation:
      `E[B prod R_i] >= 2E[B] - E[B e^(L(1-z))] >= E[B]/25 >= 1/25`.
  - *Checks.*
    - Each step was re-derived by hand.
    - The boundary bounds hold on a 200,001-point grid for labels `<= 30`
      (worst margin x1.33).
    - The cutoff conclusion holds exactly, and the agent verifier reruns.
  - Example: the cutoff for `(3,3,2)` drops from 56,000 (LR) to 560.
- **Remaining region (FM-MECH13).**  After THREE, LS, LS applied after the
  parity transformation, LR and LR4, the uncovered H-only consumer
  instances satisfy three strict inequalities (note entry of FM-MECH13).
  The region is still infinite.  For example `phi_n(h_(2n+1)^2 h_(2n)
  h_1^(2n))`, `n >= 4`, lies outside all criteria, while its exact values
  are positive for `n = 4..7`.
- **Theorem BW (a balanced-wedge slice; FM-MECH14, astra_max_ceres;
  verified by the main agent, `mech/wedge_repro.py`).**  For `m >= r >= 2`:
  `phi_r(h_(2m+1)^2 h_(2m) h_1^(2r)) > 0` and `phi_r(h_(2m)^3 h_1^(2r)) > 0`.
  The other parity patterns of three labels from `{2m, 2m+1}` vanish.
  - *Proof.*
    - The word splits exactly into a one-label part `A >= 2 tau_1`, which
      is nonnegative by Theorem OL, minus corrections `C(K)`.
      Telescoping `W(p,q) = F_p(q) - F_p(q+2)` over the full fusion
      intervals gives these corrections, with `K = 2r+m`.
    - The corrections are bounded by `24 T_r^2` or `27 T_r^2`, a tail
      binomial.
    - `b_r / T_r > 3r/4` then gives
      `2 tau_1 = 8(2r-1) b_r^2/r^2 > 27 T_r^2`.
    - Definitions (FM-CHK27).
      - `c_j = [z^j](1+z)^3(1-z)^(2r-3)` and `N = 4r-3`.
      - `W(p,q) = B_i c_j - c_i B_j`, with `i = (N+p+q)/2+1` and
        `j = (N+p-q)/2`.
      - `F_p(q) = c_(i-1) c_j - c_i c_(j+1)`, `b_r = C(2r-3, r-1)` and
        `T_r = C(2r-3, floor((3r+1)/2) - 1)`.
      - The corrections are
        `C_oo = c_K^2 - c_(K-1)^2 + 2 c_K (c_(K+1) - c_(K-1))` and
        `C_ee = 3(c_K^2 - c_(K-1)^2)`, with `K = 2r+m`.
    - Small levels, corrections as (two-odd, all-even):
      - `r = 2`: `m = 2` gives `(-1, -3)`; `m >= 3` gives `(0, 0)`.
      - `r = 3`: `m = 3` gives `(-14, -24)`, `m = 4` gives `(-1, -3)`;
        `m >= 5` gives `(0, 0)`.
      - `r = 5`: `m = 5..8` give `(-568, -1152)`, `(23, -21)`,
        `(-14, -24)`, `(-1, -3)`; `m >= 9` gives `(0, 0)`.
      - All are within the caps.
  - The agent verifier reruns.  The main agent recomputed
    `phi_4(h_9^2 h_8 h_1^8) = 4440` and
    `phi_5(h_11^2 h_10 h_1^10) = 60434` directly.
  - Scope: a two-parameter slice (specific labels, `a = 2r`), so a family
    result under the full-cone rule.  It covers FM-MECH13's uncovered
    family (`m = r = n`).
  - Next comparison, for general labels at `a = 2r`:
    `tau(U_A U_B U_C) >= W(U_A U_B, U_C) + W(U_A U_C, U_B) + W(U_B U_C, U_A)`.
- **Theorem LL (three factors with large labels, every level and
  suffix; FM-MECH15, astra_max_ceres; verified by the main agent,
  `mech/ll_repro.py`, `mech/ll_check.py`).**
  - If `min(u,v,w) >= a + 2r - 4`, or (sorted `u >= v >= w`)
    `u - v + w >= a + 2r - 3`, then `phi_r(h_u h_v h_w h_1^a) >= 0`.
  - *Proof.*
    - With `e = 2r-3` (odd), `N = a+e`, `A,B,C = u+1, v+1, w+1`, the
      divided-difference expansion gives
      `phi = tau(U_A U_B U_C) - W(U_A U_B, U_C) - W(U_A U_C, U_B) -
      W(U_B U_C, U_A)`.
    - `tau(U_l) = phi_(r-1)(h_(l-1) h_1^a) >= 0` by Theorem OL, and the
      Clebsch--Gordan multiplicities are `>= 0`.
    - `W(p,q) = 0` for `p + q > N`.  Expand the weight, which is
      homogeneous of degree `N`, in monomials `x^i y^(N-i)`.  `U_p(x)`
      is orthogonal to `x^i` for `i < p`, and `U_q(y)` is orthogonal to
      `y^(N-i)` for `N - i < q`.  So every monomial pairs to zero when
      `p + q > N` (explicit form requested by FM-CHK28).  At the boundary, `W(p,N) = -[p=0]`, so each surviving
      correction contributes `+1` when two labels match.
  - *Level 1* (FM-CHK28 repair): `phi_1(h_u h_v h_w h_1^a)` is an
    invariant multiplicity, so it is `>= 0`; the proof above is for
    `r >= 2`.
  - *Checks.*  The main agent re-derived the expansion and the support
    lemma, reran the agent verifier, and checked the exact formula in 198
    random large-label cases with its own kernel code, all agreeing.
  - FM-CHK28 (luna_max_venus, own code): ACCEPT after these two repairs.
    It ran 1,411 divided-difference identities and 492 region cases, and
    confirmed that the leftover region is the exact complement.  It notes
    that the case `w = 0` is a two-label word, where the strip and
    Theorem RF of `(E)` apply.
  - Left over, among three-factor words: `w <= a+2r-5`,
    `u - v + w <= a+2r-4`, and below LR4's cutoff.  This region is still
    infinite, since the two largest labels can grow together.
- **Three-factor closed form and Theorem T3R (FM-MECH16,
  astra_max_ceres; verified by the main agent, `mech/closed3_repro.py`,
  `mech/closed3_check.py`).**
  - *Closed form, all labels.*  Sort `A >= B >= C` (`= u+1, v+1, w+1`), and
    put `e = 2r-3`, `N = a+e` and `alpha = (N+A-B-C)/2`,
    `beta = alpha+C`, `gamma = alpha+B`, `delta = alpha+B+C`.  For
    `r >= 2` and `a + u + v + w` even (odd parity gives 0; `r = 1` is an
    invariant multiplicity, so it is `>= 0`; FM-CHK29 scope repair),
    `phi_r(h_u h_v h_w h_1^a) = M - R`, where:
    - `M = sum_(k=alpha)^beta D_k - sum_(k=gamma+1)^(delta+1) D_k >= 0`
      (by Theorem OL, via the Clebsch--Gordan expansion of
      `tau(U_A U_B U_C)`);
    - `R = Q(c_alpha, c_beta, c_gamma, c_(delta+2))
      - Q(c_(alpha-1), c_(beta+1), c_(gamma+1), c_(delta+1))`, with
      `Q(x,y,z,t) = (x+t)(y+z) - xt - yz`.
    - The three cross terms telescope via
      `W(U_p, U_q) = F_p(q) - F_p(q+2)`.
  - *Checks.*  Direct values match in 1,176 cases, and main/cross
    separately in 10,890 (agent code, rerun).  The main agent's kernel
    evaluator agrees in 297 random cases.
  - **Theorem T3R (three new regions, every level).**
    - (1) `a = 2r-3` with exactly one odd label.  The weight
      `(x^2-y^2)^e` is even, so both remaining corrections vanish, and
      each mixed two-label word is `>= 0` by the strip of `(E)`.
    - (2) `a = 2r-3`, all labels odd, and the two smaller `== 1 mod 4`.
      Ordered binomial magnitudes give `R <= 0`.
    - (3) `u + v - w >= a + 2r - 3` together with `|a - (2r-3)| <= 1` or
      `min(a, 2r-3) <= 2`.  The corrections vanish by support, and the
      terms are proved cases of `(E)`.
  - *Obstruction to an abstract proof.*  The artificial antisymmetric
    sequence `(1,0,-1,0,1,0,-1,0,1,0,-1)` satisfies OL and every `(E)`
    inequality, yet gives `M - R = -1` at `A = B = C = 4`.  So the band
    needs structure specific to `(1+z)^a (1-z)^e`.  The knob violated is
    abstractness, which the consumer does not need.
  - Open: `M >= R` for the actual binomial coefficients on the rest of the
    band.  The first residual is `a = e` with labels `== 3 mod 4`.
  - FM-CHK29 (luna_max_venus, own code): ACCEPT after the scope repair.
    - Identities checked: 1,764 closed-form/direct, 14,520 endpoint
      sums, 7,645 telescoping, 3,360 separate main/cross.
    - T3R regions checked: 1,176, 1,792 and 6,839 cases.
    - The artificial-sequence obstruction is confirmed.
- **Theorem DS (the full diagonal strip; FM-MECH17, astra_max_ceres;
  verified by the main agent, `mech/ds_repro.py`, `mech/ds_check.py`;
  independent check FM-CHK30 (luna_max_venus): ACCEPT).**
  - `phi_r(h_u h_v h_w h_1^(2r-3)) >= 0` for every `r >= 2` and all
    `u, v, w >= 1`.
  - Also, for odd `u, v` and even `w`:
    `phi_r(h_u h_v h_w h_1^(2r-2)) = phi_r(h_u h_v h_(w-1) h_1^(2r-3)) +
    phi_r(h_u h_v h_(w+1) h_1^(2r-3)) >= 0`.
  - *Proof, all-odd case* (the other parities are T3R(1) or zero).
    - At `a = e`, `c_(2h) = (-1)^h C(e,h)` and the odd coefficients vanish.
    - Ordered binomial magnitudes `x >= y >= z >= t` give `R <= Rbar`,
      where `Rbar = I + J -+ K`, `I = xy - zt`, `J = xz - yt`,
      `K = yz - xt`.
    - Tail inequality: `b_(k-2) + b_(k+1) >= 2 b_k` for
      `k >= (e+5)/2`.
    - A scalar comparison lemma settles the base case `S = 1` (smallest
      half-label 1).  The exceptional pair `(P,Q) = (2,1)` has the closed
      form `C(2n+1,n)^2 * 8(n+1)(n^2+3n+9)/[(n+2)^2 (n+3)^2]`.
    - Increment lemma: raising the two smaller half-labels together
      cannot decrease `M - Rbar`.  The increment is factored exactly and
      is bounded in two ratio ranges, with worst case `31/54`.
  - *Checks.*
    - Agent verifier (rerun): `R <= Rbar <= M` in 15,960 cases;
      increment identities in 13,566; base identities in 2,394;
      Catalan bridge in 420 cases.
    - Main-agent kernel grid: 14,690 values on `a = 2r-3`
      (`r <= 11`, labels `<= 25`) and 5,808 on the adjacent branch, none
      negative.
  - **Theorem AS (both adjacent strips; FM-MECH18, astra_max_ceres;
    verified by the main agent, `mech/as_repro.py`, `mech/as_check.py`).**
    - `phi_r(h_u h_v h_w h_1^a) >= 0` for `a = 2r-4` and `a = 2r-2`, every
      `r >= 2`, all labels.  With DS, every three-factor word with
      `|a - (2r-3)| <= 1` is proved.
    - *Ingredients.*
      - A binomial ratio monotonicity lemma, via the `tanh`/`atanh`
        composition `F -> (F+c)/(1+c k^2 F)`.
      - Pascal-identity forms of `I, J, K`, giving `I >= K >= 0` and
        `J >= K`.
      - A joint increment for all-even labels, a polynomial in the gaps
        `q_i >= 0` whose 20 coefficients are all `>= 2`.
      - The DS bound extended to both row parities, with three exceptional
        base gaps in closed form.
      - A transfer for two odd labels on the lower strip:
        `D = D_0(k) + D_0(k-1)`, `M = M_- + M_+`, `R <= U_- + U_+`.
    - *Checks.*
      - The agent verifier reruns: 14,000 bounds, 5,691 increments, 2,548
        transfers and 1,200 Catalan bridges.
      - The main-agent kernel grid on both strips (`r <= 10`, labels
        `<= 22`) gave 9,108 + 9,108 values, none negative.
    - FM-CHK30 (luna_max_venus, own code): ACCEPT.
      - Checks run: 87,125 ratio values, 4,608 + 4,608 Catalan/endpoint
        cases, 2,860 DS bounds, 4,368 all-even cases, and 2,548 + 2,548
        transfers.
      - It also wrote out the `rho >= 3/2` step of the DS increment.
    - Next: `a = e+2`, with `c_(2h) = (-1)^h [C(e,h) - C(e,h-1)]` and
      `c_(2h+1) = 2(-1)^h C(e,h)`.  The comparison and the increment need
      replacing when consecutive magnitudes differ.
  - **Offset induction: the whole three-factor sector reduces to one
    inequality (P) (FM-MECH19, astra_max_ceres; reproducer rerun by the
    main agent, `mech/pinduct_repro.py`).**
    - *Fold.*  Under `a <-> e`, `M` is invariant and `(R_C, R_B, R_A)`
      picks up the signs `((-1)^C, (-1)^B, (-1)^(B+C))`.
    - *Pascal increments.*  For `a >= e`, the quantities
      `D_k = c_k^2 - c_(k-1) c_(k+1)`, `E_k = c_k c_(k-1) - c_(k-2) c_(k+1)` and
      `H_k = c_k^2 - c_(k-2) c_(k+2)` are all `>= 0`.  This uses Newton's
      inequalities for the real-rooted even and odd subsequences
      `(1-w)^e sum C(d, 2j (+1)) w^j`.
    - Under `q_k = c_(k-1) + 2c_k + c_(k+1)`, the effect of `a -> a+2`,
      `D[q]_k - D[c]_k = L_k`, with
      `L_k = D_(k-1) + D_k + D_(k+1) + 2E_k + 2E_(k+1) + H_k >= 0`.
    - *Reduction.*
      - With the consumer margin `G_(e,a) = M - R_(e,a)`, taken over the
        physical orientation(s):
        `G_(e,a+2) - G_(e,a) = sum_(k=alpha)^beta L_k -
        sum_(k=gamma+1)^(delta+1) L_k - [R(q) - R(c)]`.
      - So **(P)**, `R(q) - R(c) <= sum_(alpha..beta) L - sum_(gamma+1..delta+1) L`,
        together with the base strips DS and AS, proves every three-factor
        word at every level by induction in `a` (in steps of 2 on each
        side).
    - *Checks.*  (P) holds at all 40,950 tested steps (`e <= 16`,
      `a <= e+17`, labels `<= 12`; agent code, rerun), with 1,488
      independent orientation bridges.
    - KILL (single-label Pascal step): `phi_3(h_4^2 h_6 h_1^4) = 60`, below
      each neighbouring sum `62, 62, 66`.  The knob violated is a stronger
      induction, which the consumer does not need.
    - FM-CHK31 (luna_max_venus, own code): ACCEPT for the reduction
      (fold, Newton/real-rootedness, Pascal identities, orientation
      bookkeeping, base cases).
      - Its independent (P) screen covered 250,250 triples (`e <= 10`,
        `a <= e+40`, labels `<= 20`), with no failure; the minimum slack
        is 0, so equality cases occur.
      - It also ran 39,325 identity checks.
    - FM-SEC10 (luna_max_venus): the support edge of (P).
      - Zero slack occurs exactly when `alpha >= N+2` (both windows and
        both corrections vanish).  Slack 1 occurs exactly when
        `alpha = N+1`.
      - (P) is proved for every `alpha >= N-1`; at `alpha = N-1` the
        slack is an explicit positive polynomial in `N` and
        `d = a - e`.
      - In the interior `alpha <= N-2` the screened minimum slack is 7
        (250,250 triples).
    - FM-MECH20 (astra_max_ceres): no proof of (P); an exact form and a
      stronger obstruction.
      - *Square balance.*
        - For each orientation `sigma`,
          `2F_sigma = E(c) + X_sigma^2 - Y_sigma^2`.
        - Here `E = sum_(alpha..beta) g_k^2 - sum_(gamma+1..delta+1) g_k^2`,
          and `g_k = c_(k+1) - c_(k-1)` are the coefficients of
          `(1+z)^(a+1)(1-z)^(e+1)`.
        - `X_sigma`, `Y_sigma` are signed four-endpoint sums, and they
          factor through shift operators.
      - All three correction increments are telescoped sums of one
        determinant kernel `J(i,j) = v_i q_j - q_i v_j`, with
        `v = (T+1)^2 c`.
      - KILL (abstract induction, strong form).
        - The sequence `(3,0,-5,0,5,0,-5,0,5,0,-3)` satisfies all of:
          `D, E, H >= 0` for `c` and `q`; right-half monotonicity; every
          old `(E)`; every old three-factor margin `>= 0`.
        - Yet the margin goes from 1 to -15 at `(3,3,3)`.
        - The knob violated is abstractness.  A proof must use the
          coefficient recurrence
          `(k+1) c_(k+1) = (a-e) c_k - (N-k+1) c_(k-1)`, or the new-offset
          `(E)`; the artificial `q` fails the latter.
    - FM-SEC11 (luna_max_venus): (P) proved also when `alpha <= N-2` and
      `beta >= N+2`; there every endpoint pair has an index `>= beta`,
      so the corrections vanish by support.
      - Window identity: the right side equals `sum_k w_k (L_k - L_(k+1))`
        with explicit trapezoid weights `w_k`.
      - Tightness: in 245,379 screened interior cases the correction
        reaches 97.4% of the right side (minimum relative slack 2.6%,
        at `(e,a;u,v,w) = (0,39;6,1,1)`), so a crude bound cannot
        work.
      - KILL (right-half log-concavity).  At `e = 2`, `a = 5`, the
        coefficients `[1,3,1,-5,-5,1,3,1]` have non-unimodal right-half
        magnitudes, although (P) holds there with slack 164.
    - FM-SEC12 (luna_max_venus): symmetric unimodality of
      `L_k = D[q]_k - D[c]_k`.
      - Proved for `e = 0`, where the difference factors with a
        positive quartic `P(s,k)`, and for `a = e`, where
        `L_(2j+1) - L_(2j) = (b_(j+1) - b_j)(b_(j-1) + b_j + b_(j+1))`.
      - Screened on 50,601 more rows, with no failure.
      - Unimodality would make the right side of (P) nonnegative by
        window ordering.
      - The square-balance components have mixed signs in the residual
        region, so termwise positivity fails there.
    - Open: (P) for `alpha <= N-2`, `beta <= N+1`.
  - **Theorem G3 (interior regions, directly, without (P); FM-MECH21,
    astra_max_ceres; verified by the main agent, `mech/g3_repro.py`,
    `mech/g3_check.py`).**
    - Assume (G0): `C` is even, `0 <= alpha < beta <= N` and `gamma >= N+1`.
      In the three-factor closed form `phi = M - R`, the coefficients with
      indices `gamma, gamma+1, delta+1, delta+2` vanish (`delta = gamma+C`,
      `C >= 2`).  Hence the upper `D`-window is zero and
      `R = c_alpha c_beta - c_(alpha-1) c_(beta+1)`; this endpoint product
      stays in the formula (repair from FM-CHK32).  With `C = 2h` and
      `g_k = c_(k+1) - c_(k-1)`:
      `phi_r(h_u h_v h_w h_1^a) = sum_(1<=i<=j<=h) [g_(alpha+2i-1) g_(alpha+2j-1)
      - g_(alpha+2i-2) g_(alpha+2j)]`.
    - `g` is the coefficient row of the real-rooted
      `(1+z)^(a+1)(1-z)^(e+1)`.  By Newton, on any block where `g` keeps one
      sign every minor is `>= 0`, and the `i = j = 1` minor is `> 0`.  The
      fold `a <-> e` preserves each summand.
    - One-sign blocks come from `K_(m+1)(X; N+2)`, `m = min(a,e)`, having no
      root in `(d, Y)`, where `d = A-B-C` and `Y = A-B+C`.  This holds in:
      - the central region, `m` odd and `(m+1) Y^2 <= 2(N-m+3)`
        (reciprocal-root bound);
      - the outer region, `d >= 0` and `d^2 >= 4m(N-m+3)` (Jacobi-matrix
        bound).
    - It covers families outside all earlier cutoffs, e.g.
      `phi_3(h_18^2 h_3 h_1^29) = 64,110,062,011,500`.
    - *Checks.*  The agent verifier reruns.  The main agent's kernel
      evaluator confirmed the minor identity in 300 random (G0) cases, and
      positivity in the 16 of them inside `(GC)` or `(GO)`.
    - Open: (G0) windows crossing a root of `K_(m+1)`; `C` odd or
      `gamma <= N`.
  - **Theorems G3X and G3O (root crossing, odd `C`, part of `gamma = N`;
    FM-MECH22, astra_max_ceres; verified by the main agent,
    `fm39/g3x_repro.py`, `fm39/g3x_check.py`, `fm39/g3x_region_check.py`).**
    - *Decomposition, all three corrections present.*  Put
      `A_x = c_x - c_(x+C)`, `B_x = c_(x-1) - c_(x+C+1)` and
      `P_C(x) = sum_(k=x..x+C) D_k - c_x c_(x+C) + c_(x-1) c_(x+C+1)`.
      Then, with `l = alpha` and `tau = gamma+1`,
      `phi_r(h_u h_v h_w h_1^a) = P_C(l) - P_C(tau) - (A_l B_tau - B_l A_tau)`.
      For `C = 2h`, `P_C(x)` is the minor sum `T_C(x)` of G3.  For every
      `C` there is a boundary-free form in `f_k = c_k - c_(k-1)` (the row
      of `(1+z)^a (1-z)^(e+1)`):
      `P_C(x) = sum_(s=2..2C) [f_(x+floor(s/2)) f_(x+ceil(s/2))
      - f_(x+max(0,s-C-1)) f_(x+min(s,C+1))]`.
    - *Single-root lemma.*  If `F = (x - lambda) H` with `H > 0` and
      `-kappa <= (log H)'' <= 0` on `[L, U]`, and
      `kappa (lambda - L)(U - lambda) <= 1`, then every balanced product
      dominates its outer product on `[L, U]`, across the root.  The proof
      has three cases (repair from FM-CHK33):
      - the outer pair does not cross `lambda`: concavity of `log|F|` on
        that side;
      - the outer pair crosses and the inner pair does not: signs;
      - both pairs cross: the curvature estimate
        `H(x+s) H(y-s) Delta [1 - kappa (lambda - x)(y - lambda)] >= 0`,
        with `Delta = s(y-x-s)`.
      The Gamma-interpolation bound `0 < -(log B_n)'' <= n/(n^2 - X^2)`
      needs `|X| < n`.  G3X uses it only for `Y^2 <= n/6`.  The
      reciprocal-root bound `sum xi^-2 <= s/(n-2s+2)` for the positive roots
      of `K_(2s)` needs `n >= 2s`.
    - **Theorem G3X.**  `C` even, `m = min(a,e)` even and `>= 2`,
      `0 <= alpha < beta <= N`, `gamma >= N+1` and `3 m Y^2 <= N - m + 4`
      (`Y = A - B + C`).  Then `phi > 0`; the window may cross the
      central root of `K_(m+1)`.  Example: `(r,a;u,v,w) = (150,4;153,152,3)`.
    - **Theorem G3O.**  `C` odd, `0 <= alpha < beta <= N`, `gamma >= N+1`
      and `(e+1)(Y+1)^2 <= 2(a+2)`.  Then `phi > 0`.  Example:
      `phi_3(h_18^2 h_2 h_1^30) = 147,831,720,921,000`.
    - *`gamma = N`.*  There `phi = P_C(alpha) + c_alpha - c_beta`
      (`c_N = -1`, `A_(N+1) = 0`, `B_(N+1) = -1`, `P_C(N+1) = 0`).  Both
      theorems extend, with their other hypotheses kept: G3O when `r` is
      even, G3X when `alpha + m/2` is even.  In the G3X branch `a = m`, so
      `N` is odd; with `C` even and `alpha` integral, `A - B` is odd, hence
      `A - B > 0`, which the sign argument uses.
    - *Checks.*  The agent verifier reruns.  The main agent's direct split
      evaluator (definition-level ballot moments, cross-checked against
      `fm3kern.py`) confirmed the decomposition in 400 random cases and
      positivity in 9,513 (GX) words (every one with `d < 0`, so the window
      crosses `X = 0`), 21,316 (GOdd) words and 282 `gamma = N` words.
    - FM-CHK33 (luna_max_venus, own exact code):
      - ACCEPT: G3O, the reciprocal-root lemma, and the example families
        and values.
      - REPAIR (applied above): the Lemma 2 case split; the domains
        `|X| < n` and `n >= 2s`; and the `A - B > 0` step in the
        `gamma = N` G3X branch.
      - Presentation repairs to the agent text: the square-balance
        citation, and the expansion from `phi = M - R` to the
        decomposition.
      - Remaining after the repairs: `P_C(alpha) - P_C(gamma+1) >=
        A_alpha B_(gamma+1) - B_alpha A_(gamma+1)` outside these regions.
    - Kill (knob: termwise sign of the minors, no consumer use): single
      minors can be negative, e.g. `a = e = 3`, `(5,5,3)` gives minors
      `24, -16, 24`.
  - **The whole G0 branch reduces to one real-rootedness inequality
    (main agent, `fm39/g0_scan.py`, `fm39/window_explore.py`,
    `fm39/window_realrooted.py`, `fm39/window_f.py`).**
    - On `gamma >= N+1`, `phi = P_C(alpha)`.  All 533,190 G0 words with
      `r <= 8`, `a < 40` (both parities of `C`) are positive.
    - Conjecture W: for every real polynomial `P = sum c_k z^k` with only
      real roots (any signs) and every window, `P_C(x) >= 0`.  Screens:
      523,745 windows of `(1+z)^a (1-z)^e`, `a, e <= 30`; 1,500 random
      real-rooted rows (all windows); 95,659 windows in the `f`-form
      `S(w) = sum_(j=1..M-1) w_j^2 + sum_(j=1..M-2) w_j w_(j+1)
      - w_0 sum_(s=2..M) w_s - w_M sum_(s=1..M-2) w_s`, for arbitrary
      real-rooted rows; no failure.
    - Controls: one complex root pair gives negative windows in 138 of
      300 rows; random integer rows in 197 of 300.  The `C`-increment
      `P_(C+1) - P_C` and the shifted sums `sum_k T(k,k+s) - T(x,x+C+s)`
      fail often.
    - Conjecture W implies the entire G0 branch of the three-factor
      stratum at every level.
  - **The `gamma <= N` branch (main agent, `fm39/gammaN_branch_scan.py`,
    `fm39/mirror_check.py`).**
    - *Mirror form.*  For anti-reciprocal rows, `P_C(tau) = P_C(l - A - 1)`
      and
      `phi = P_C(l) - P_C(l-A-1) - (g_l + g_(tau-A-1)) ^ (g_tau + g_(l-A-1))`,
      with `l = alpha`, `tau = gamma+1`.  Checked in 730 cases.
    - *Size and margin.*  Outside the strips (`|a-e| >= 2`, `a >= 1`),
      `r <= 6`, `a < 30`: 50,370 words, all positive.  Throughout,
      `phi / tau(U_A U_B U_C) >= 0.46`.  The tightest words are
      `h_2^3 h_1^a` at `r = 2`.
    - At fixed labels the ratio decays like `1/a^2`, e.g.
      `0.073, 0.019, 0.0048` at `a = 80, 160, 320` for `h_2^3`, `r = 2`.
      That range is inside Theorem LS.
    - Assigned as FM-SEC14 (luna_max_mars).
    - *FM-SEC18 (luna_max_uranus): low-`s` layers.*  With `p = A - C`,
      `q = B - C` and `s = N - gamma`, one has `N = C + p + q + 2s` and
      `alpha = p + s`.
      - For `s <= 2`, exactly `phi = P_C(p+s) + R_s`, with `d = a - e`:
        - `R_0 = c_p + c_q`;
        - `R_1 = d(c_(p+1) + c_(q+1)) - (c_p + c_q) - 1`;
        - `R_2 = c_2 (c_(p+2) + c_(q+2)) - d(c_(p+1) + c_(q+1)) - 1 - (d^2+N)/2`,
          where `c_2 = (d^2 - N)/2`.
      - Conditional slice: under the listed (E) instances,
        `phi_r(h_(C-1)^3 h_1^(C-e)) >= 2` for `C >= max(3, e)`.
      - `R_s` can be negative, e.g. `(e,a,C,p,q) = (3,1,3,1,0)` has
        `R_0 = -1`, `P = 11`, `phi = 10`.  So the branch needs a
        quantitative lower bound on the window, not only W.
    - *FM-SEC30 (luna_max_uranus).*  Proved uniformly on the
      equal-label low-`s` regions, including every `s = 0`, `p = q = 0`
      case, without (E).  Bounding the combined cross wedge by the
      determinant of `H` loses too much in general, e.g. at
      `(e,a,C,p,q,s) = (1,23,3,1,0,10)`, where `phi > 0`.
    - *FM-SEC45 (luna_max_uranus): one `Q_t` wedge for the cross term.*
      - Notation: `N = C + p + q + 2s`, `U = g_(p+s) + g_(q+s)`,
        `V = g_(N-s+1) + g_(s-C-1)`, `phi = delta - U ^ V`, with
        `delta = P_C(p+s) - P_C(s-C-1) = sum_d (E + eps W)` over the
        CG labels `d` (each step an (E) value).
      - Certificate: bound each step below by `lambda_d` (`E` when
        `eps W >= 0`, else `E - ceil(sqrt R)` from `(M_i)`/`(M_j)`), put
        `Lambda = sum lambda_d`, and bound `|U ^ V|` by FM-MECH26's
        metric: `16t(4V-t)|U ^ V|^2 <= F_U(t) F_V(t)`, with
        `F_Y(t) = 4(Y_1^2 - Y_2^2) t + 4(d Y_1 - sigma Y_2)^2`.  If
        `K = 16 Lambda^2 + alpha_U alpha_V > 0`,
        `0 < L < 8VK` and `L^2 >= 4KZ` (explicit in the row), then
        `phi >= 0`.  No (E) instance is assumed.  The boundary cases
        rerun exactly.
      - Main-agent census (`fm39/gammaN_cert_census.py`): `e <= 11`,
        `a <= 20`: 33,308 of 33,480 words certified.  The failures have
        `e = 1` (`r = 2`, covered by T1) or `a <= 2`.
      - Direct form (`fm39/gammaN_metric_census.py`, exact).  At the row
        `delta` is known, so test
        `(HT3)`: some `t in (0, 4V)` has
        `16t(4V-t) delta^2 >= F_U(t) F_V(t)`.
        Over `e` odd in `3..25`, `a <= 40`, all `C >= 3`, `p >= q >= 0`,
        `s >= 0`: 1,196,247 of 1,196,426 words satisfy (HT3), 124 more
        have `U ^ V <= 0`, and 55 fail, all with `phi > 0`.  Correction
        (FM-SEC51, luna_max_uranus): the failures fall into three groups.
        - 7 with `e = 3`, `a = 37..40`, `C = 3`, e.g.
          `(e,a,C,p,q,s) = (3,37,3,15,0,11)`.
        - 17 with `a <= 2`.
        - 31 with `a >= 3` and `e >= 17`, mostly `C = 3, 4`, `p = q = 0`
          (words like `h_2^3 h_1^a` at high level), e.g.
          `(17,4,3,0,0,9)`, which is `r = 10`, `h_2^3 h_1^4`.  An earlier
          version of this entry omitted this group.
        FM-SEC51's geometric condition covers exactly the `e = 3` group:
        both W-window sweeps `<= 2 pi`, and short `psi`-arcs for the
        adverse (E)-steps.  So the split is not yet exhaustive, and there
        is also a regime `e >> a`.  With HT3, cross `<= 0` and the geometric
        condition, the `e = 3, 5`, `C = 3, 4` rows up to `a = 400` are
        covered (about 10.9M words).
      - FM-SEC67 (luna_max_uranus): sharp radial theorems
        (`fm39/sec67_radial_repro.py`, rerun by the main agent after one
        sympy API patch, `coeffs()` for a multivariate polynomial).
        - Exact joint moments: with `s = x+y`, `d = x-y`,
          `E[s^(2m) d^(2k)] = 2 (2m)! (2m+1)! (2k)! (2k+1)! /
          ((m!)^2 (k!)^2 (m+k+1)! (m+k+2)!)`.  So `phi_R(w h_1^e)` is an
          explicit hypergeometric sum for any fixed core.
        - Cutoffs (positive-coefficient numerators in
          `xi = n - R - const`), for odd `e`: `phi_R(hat S_3^3 h_1^e) > 0`
          for `e >= 2R+3`, and `phi_R(h_3^3 h_1^e)`,
          `phi_R(h_5 hat S_3^2 h_1^e)` `> 0` for `e >= 2R+1`.  For even `e`
          the value is 0 by parity (FM-CHK43 repair, e.g.
          `phi_1(hat S_3^3 h_1^6) = 0`).  LS needed `28R+41`,
          `(116R+169)/5` and `36R+53`.
        - Via the parity transform, 18 of the 31 `e >> a` failures are
          covered.  The other 13 (`h_3^3` at `R = 17..21`) and the `a = 2`
          misses are exact positive values.  So the whole screened
          `gamma <= N` box (`e = 3..25` odd, `a <= 40`; 1,196,426 words) is
          closed.
        - Open: the short-suffix side `e < 2R + O(1)` and general cores.
          For `hat S_5^3` the numerator coefficient `[xi^11] = -108(2R-9)`
          turns negative for `R >= 5`.  FM-SEC75 seeks a sharp LS for
          every word.
        - FM-SEC75/80 (luna_max_uranus): for general cores the explicit
          affine cutoff is still Theorem LS's
          `e >= (2 Lambda(w) - 2) R + 3 Lambda(w) - 4`.  At
          `e = 2R + C(w)` the concentration point sits on an edge saddle
          path, not at `s = 4, d = 0`.  The exact ratio formula closes the
          whole suffix range for `h_3^3`.  At `q = R - n = 1` the numerator
          is a positive sextic in `n`; for `q >= 2` it has positive
          coefficients in `q - 2`.
        At `e = 1` (level 2), `a <= 60`: 146,474 of 148,800 words
        satisfy (HT3), 767 have `U ^ V <= 0`, and 1,559 fail (all
        `phi > 0`, e.g. `(1,24,3,4,0,9)`).  Assigned as FM-SEC57.
      FM-SEC57 (luna_max_mars): exact closed forms on the slice `C = 3`,
      `q = 0`, i.e. `phi_2(h_(p+2) h_2^2 h_1^(p+2s+2))`, via the explicit
      row `c_k = C(a,k) - C(a,k-1)`:
      - `s = 0`: `(p^6+3p^5+7p^4+9p^3+280p^2+996p+1296)/144`;
      - `s = 1`: `(p+3)(p+4)(p^6+9p^5+35p^4-37p^3+348p^2+4972p+13680)/2880`;
      - `s = 2`: `(p+3)(p+4)(p+5)(p^7+23p^6+217p^5+725p^4-866p^3+7892p^2
        +175608p+580320)/86400`.
      Each is positive for all `p >= 0`, and the formulas are checked
      against direct evaluation for `p <= 40`.  Open: `s >= 3` (the part
      outside LS is `10s < 7p^2+51p+202`), `q > 0`, `C >= 4`.  FM-SEC65
      continues.
      FM-SEC65 (luna_max_mars): for `s >= 5` the same family factors
      exactly (`fm39/e1_q0_closedform.py`) as
      `phi_2 = B^2 (p+3)(p+4)(p+5)(p+s+8)(p+2s+4) Q(p,s) / (s^2 (s-4)^2
      (s-3)^2 (s-2)^2 (s-1)^2 (s+1)^2 (s+2)^2 (s+3)^2 (s+4) (p+2s+3))`,
      with `B = C(p+2s+3, s-5)` and `Q` an integer polynomial of degree 12
      in `p` and 10 in `s`.  So the family reduces to `Q(p,s) >= 0` for
      `p >= 0`, `s >= 5`.  `s = 3, 4` are checked directly (positive).
      Main agent (`fm39/e1_q0_Qpositivity.py`): `Q(p, 5+S)` has 15 negative
      coefficients, all at `p`-degree 3 and 4.  Its top homogeneous part is
      `p^4 (10S^4 + 20S^3 p + 16S^2 p^2 + 6S p^3 + p^4)^2`, which vanishes on
      the `S`-axis, so Polya multiplication cannot certify it.  On the
      real grid `[0,60]^2` its minimum is `3.88e12`, at the corner.  An
      AM-GM certificate with pairs alone is infeasible, and so is one with
      triangles (82,250 candidates).  Numerically `Q` is comfortably
      positive: `Q / sum |coefficient * monomial| >= 0.1355` on a log grid
      `p, S in [0, 1e4]` (minimum near `p ~ S^0.6`), with integer minimum
      `3.88e12` at `(p,s) = (0,5)`.  So a decision-procedure certificate is
      needed (discriminant / real-root isolation, or an SOS form in `p`).
      Assigned as FM-SEC71.
      **FM-SEC71 (luna_max_mars): proved** (`fm39/e1_q0_Q_certificate.py`,
      rerun by the main agent after restoring one closing parenthesis in
      the printed code).
      - Weighted-homogeneous certificate.  Group the monomials `p^i S^j`
        by weight `w = i + 2j` (the scaling `p ~ S^(1/2)`).  Then
        `Q(TY, 5+Y^2) = sum_w Y^w T^(w mod 2) H_w(T^2)`, and each `H_w` is
        of one of three kinds:
        - coefficientwise nonnegative;
        - a positive quadratic (negative discriminant) plus a
          nonnegative remainder;
        - a bounded-margin case, checked exactly.
      - The top class gives `100 S^8 (p^2 - 6S)^2 + 2400 S^10`, so
        `Q(p,5+S) >= 2400 S^10 > 0`.  Also `Q(p,5)` has positive
        coefficients.
      - With the exact small-`s` forms, `phi_2(h_(p+2) h_2^2 h_1^(p+2s+2))
        > 0` for all `p, s >= 0`.  This closes the `C = 3`, `q = 0` family
        of the `e = 1` branch.
      - Screen: `Phi(C,p,q,s) >= Phi(3,p,0,s)` on 22,680 cases
        (`3 <= C <= 12`, `p <= 15`, `q <= 8`, `s <= 20`), minimum 0
        (`fm39/e1_comparison_screen.py`).  If proved, it closes the whole
        `e = 1` branch.  FM-SEC81.
      - FM-SEC81 (luna_max_mars) reduces the comparison to two step
        inequalities, `Phi(C,p,q+1,s) - Phi(C,p,q,s) >= 0` and
        `Phi(C+1,p,q,s) - Phi(C,p,q,s) >= 0`.  They are screened on 57,750
        and 60,060 cases (minima 30 and 8) but not proved.  FM-SEC84.
      - FM-SEC84 (luna_max_mars): no proof of either step.  For
        `s >= C+2`, normalizing by `B = C(N, s-C-2)` gives
        `Phi = B^2 F(C,q;p,s)`, where `F = A_C(h) - A_C(l) - U ^ V` is an
        explicit area form in the ratios `d_j`.  The steps become rational
        brackets with exact binomial ratio factors.  No uniform positivity
        decomposition was found, and the cases `s < C+2` are separate.
        Screens rerun: 60,060 and 57,750 cases, minima 8 and 30.  The
        `e = 1` branch stays open at the two steps.
    - Strict OL, empirically (`fm39/strict_ol_scan.py`): on every row
      `a, e <= 100`, all 512,600 steps with `2k > N` have `delta_k > 0`,
      so there are no flat steps.  FM-SEC64 is proving it.
    - *FM-SEC53 (luna_max_mars): suffix `a in {1,2}`.*  The G0 branch
      there follows from Theorem T3R (`min(a,e) <= 2`), and `r = 2, 3`
      follow from T1-5 and R3-5.  Under the parity transform the rest
      maps to level-1 words with two or three `hat S` and `h_1^e`
      (FM-SEC52's target), or to three-factor words at `e = 1` with a
      long suffix (FM-SEC57).  Residual families on the `gamma <= N`
      side: `phi_((e+3)/2)(h_e h_2^2 h_1)` (one short of T3R's spread;
      63 at `e = 5`) and `phi_((e+3)/2)(h_2^3 h_1^2)` (59, 27144, 291108
      at `e = 5, 11, 13`).
    - **Binary forms for three factors (main agent,
      `fm39/phi3_binary_form.py`, `fm39/phi3_binary_census.py`).**
      - Through the recurrence, `phi_r(h_u h_v h_w h_1^a)` is a binary
        quadratic form in two consecutive coefficients.  Its coefficients
        are rational in `(a, e, labels)`.
      - Definiteness is independent of the base point, since a change of
        base point is an invertible linear change of variables.
        Definiteness implies `phi >= 0` for every solution of the
        recurrence, and in particular for the actual row.
      - Census, exact rational arithmetic, over the box `2 <= r <= 8`,
        `1 <= a <= 40`, `|a-e| >= 2`, `2 <= w <= 14`, `w <= v <= N+2` and
        `v <= u <= v+w+1`: all 373,673 words give a positive definite form
        (133,210 on the `gamma <= N` branch, 240,463 on G0).  The form
        reproduces `phi` exactly in the sanity checks.
      - Every word in that box therefore has a ratio-free certificate,
        including every `gamma <= N` word there.  On a finite box this is
        no stronger than exact evaluation of `phi`.  Its significance is
        structural: it points to definiteness as a uniform mechanism.
      - It is not universal beyond three `h` factors (FM-SEC36,
        luna_max_uranus, and the main agent's `fm39/word_binary_census.py`).
        - Three-label `hat S` patterns are definite or PSD on
          FM-SEC36's grid (`r <= 4`, `a <= 6`, labels `<= 5`).
        - Larger grid (`r <= 5`, `a <= 20`, labels `<= 9`): 36,396 words,
          of which 8 are not PSD.
          - `h h hat S`: 6 of 14,976, all at `r = 1`, e.g.
            `h_2^2 hat S_9 h_1^9`.
          - `h hat S hat S`: 1 of 15,120, namely `h_8 hat S_3^2` at `r = 5`,
            `a = 0`.
          - `hat S^3`: 1 of 6,300, namely `hat S_2 hat S_3 hat S_9` at
            `r = 5`, `a = 0`.
        - Four-factor words have indefinite forms, e.g.
          `h_2 hat S_2^3 h_1^4` at `r = 1`, yet `phi = 75 > 0` on the actual
          row.
        - So beyond three `h`'s the actual coefficient ratio must be used,
          e.g. through ratio bounds as in Theorem OL's steps 2-5.
        - For the W windows (FM-SEC37, luna_max_neptune) PSD fails on the
          `e = 3` support edge for every `a >= 21` at `C = 3` and every
          `a >= 24` at `C = 4`.  The actual values are positive.  There
          `D_C >= C^2`, so Theorem LD plus the energy bound covers the
          tails via `(sqrt D_C - 1)(sqrt D_C - C) >= 0`.
        - Even for three `h`'s it fails at large `a` (FM-SEC40,
          luna_max_uranus; `fm39/eqlabel_definite_repro.py`).  The form is
          proved positive definite on the infinite families
          `u = v = w = 2`, `a = e +- 1`.  At fixed labels it fails first
          at `r = 2`, `a = 272`, `(u,v,w) = (2,2,2)`, where `phi > 0`.
      - **Kill (FM-MECH27, astra_max_ceres; `fm39/binary_form_kill_repro.py`).**
        Uniform definiteness is false.
        - For `u = N + v + w - 1` (any `v >= w >= 2`) the form is
          `(p - dq/4)^2 + [8(N+2) - d^2] q^2/16` in `(p,q) = (y_(N-1), y_N)`.
          It is indefinite once `(a-e)^2 > 8(a+e+2)`.
        - It also fails without support truncation, e.g. at
          `(r,a,u,v,w) = (2,36,26,2,2)` using only interior indices.
        - Knob violated: positivity on every solution of the interior
          recurrence, which the consumer never uses.
        - The actual row satisfies the endpoint condition `y_1 = d y_0`,
          which selects the actual row up to scale.  The boundary-aware
          identity `F = ((d^2+N+2)/2) q^2 + (p - dq)(p + dq/2)` gives
          `phi = (d^2+N+2)/2 > 0` there.
        - These words have `u > v+w+1`, outside the census box, and lie in
          the support region of Theorem LL.  Within
          `u <= v+w+1` no failure is known.
        - So definiteness is a sufficient device on part of the cone; a
          uniform proof must keep the endpoint information.
    - *Reduction to the basic inequalities, partial (main agent,
      `fm39/basis_lp.py`, `fm39/basis_lp2.py`).*
      - As quadratic forms on anti-reciprocal sequences, the LP tested
        whether each three-factor value on this branch is a nonnegative
        combination of:
        - OL differences `D_k - D_(k+1)` (`2k > N`);
        - W windows `P_C(x)`;
        - the (E) forms `D_j - D_i +/- W_ij`;
        - and all their twists.
      - Result at `r <= 3`, `a <= 9`: 244 of 322 words are such
        combinations, and 78 are not.  The first is `h_2^3 h_1^4` at
        `r = 2`.
      - So W, (E) and OL alone do not give the whole `gamma <= N` branch
        in this linear sense.
  - **Theorem G0E: the G0 branch follows from the two-label inequality
    (E) (main agent, `fm39/w_from_e.py`, `fm39/window_E.py`).**
    - *Statement.*  Let `e = 2r-3`.  If `gamma >= N+1`, then
      `phi_r(h_u h_v h_w h_1^a) = sum_(d in CG(A,B)) [tau(U_d U_C) - W(d,C)]`.
      Each summand equals `D_j - D_i - W_ij` for the pair `(j,i)` of the
      labels `(d, C)`, i.e. the two-label value.  So (E) at odd `e` implies
      `phi >= 0` on the whole G0 branch, for every `r` and `a`.
    - *Proof.*
      - In the split identity `phi = tau(U_A U_B U_C) - W(U_A U_B, U_C) -
        W(U_A U_C, U_B) - W(U_B U_C, U_A)`, the last two cross terms
        vanish on G0 by support.
      - Expanding `U_A U_B = sum_d U_d` and using
        `tau(U_p U_q) = D_j - D_i` (telescoping of Theorem OL's one-label
        values) gives the formula.
      - For `d < C` the pair has `j < N/2 < i` with `i + j >= N+1`.  The
        mirror `D_j = D_(N-j)`, `W_ij = (-1)^e W_(i,N-j)` (from
        `c_(N-k) = (-1)^e c_k`) maps it to (E) at `(N-j, i)`, which lies in
        (E)'s range since `N - j < i`.  QED.
    - *Window form.*  Equivalently, the identity
      `P_C(x) - P_C(x+1) = D_x - D_(x+C+1) - W_(x+C+1,x)` (symbolic)
      telescopes to
      `P_C(x) = sum_(y=x..N) [D_y - D_(y+C+1) - W_(y+C+1,y)]` for
      `2x >= N - C`.  That is the consumer range of W, since `A >= B`.
      Checked on 63,350 base windows, with every step `>= 0`; the mirror
      identity was checked on 32,620 pairs.
    - *Unconditional consequences.*  The G0 branch holds wherever every
      needed pair `(d, C)`, `d in CG(A,B)`, lies in a proved (E) region.
      These include `min(a,e) <= 2` (so `a <= 2` at every `r`), the strips
      `|a-e| <= 1`, and the root-free intervals of Theorem RF.
    - FM-CHK35 (luna_max_saturn, own exact code): ACCEPT on all seven
      items, conditional on (E).
      - Consumed instances, for each `d in CG(A,B)` with `d + C <= N`:
        `D_j - D_i - W_ij` at `(p,q) = (d, C)` when `d >= C`, and
        `D_j - D_i + W_ij` at the sorted pair `(C, d)` when `d < C`.
        Both come from (E)'s absolute-value form.
      - Scope: the unconditional consequences need every required pair
        to lie in a proved (E) region.  Root-free intervals do not always
        suffice, e.g. `r = 4`, `a = 3`, `(A,B,C) = (7,5,2)`, `d = 4`, where
        `(j,i) = (5,8)` crosses the root `X = 4` of `K_5`.
    - So (E) is the master inequality for the two-label stratum and the
      G0 branch alike.  W's other routes (Theorem WS, (i)+(ii), the
      long-sweep bound) prove W without (E).
    - *`gamma <= N` branch.*  Here
      `phi = sum_(d in CG(A,B)) [T - W](d,C) - W(U_A U_C, U_B) -
      W(U_B U_C, U_A)`.  Two facts are recorded:
      - the window monotonicity `P_C(l) >= P_C(l-A-1)`, a sum of (E)
        steps, held on all 51,350 screened words;
      - so did the window-(E) bound `P_C(l) - P_C(l-A-1) >= |cross|`.
      By the basis LP, however, not every such word is a nonnegative
      combination of (E) forms.
  - **Merge move and the four-wedge correction (main agent,
    `fm39/merge_check.py`, `fm39/merge_best.py`, `fm39/cross_telescope.py`).**
    - *Merge.*
      `(U_X(x) - U_X(y))(U_Y(x) - U_Y(y)) = sum_(d in CG(X,Y)) (U_d(x) + U_d(y))
      - [U_X(x) U_Y(y) + U_Y(x) U_X(y)]`.
      So two general factors become a sum of `hat S_d` factors minus a
      mixed term.  The same move applies inside any consumer word.
    - *Three factors.*  For each choice of the singleton `Z` in
      `{A, B, C}` (with `{X, Y}` the other two labels),
      `phi = sum_(d in CG(X,Y)) [T - W](d, Z) - M_Z`, where
      `M_Z = W(U_X U_Z, U_Y) + W(U_Y U_Z, U_X)`.  Each summand is a
      two-label (E) value.
    - *Telescoping.*  With `S(p,q) = g_p ^ g_q = c_p c_(q-1) - c_(p-1) c_q`,
      `j(d) = (N+d-Y)/2` and `i(d) = (N+d+Y)/2 + 1`,
      `W(U_X U_Z, U_Y) = S(j(|X-Z|), i(|X-Z|)) - S(j(X+Z)+1, i(X+Z)+1)`.
      Checked in 1,963 cases.  So `M_Z` is a combination of at most four
      phase-plane wedges.
    - *Sign census* (`r <= 7`, `a < 30`, labels `<= 15`, 65,250 words).
      Some grouping has `M_Z <= 0` in 56,458 words (87%), and then (E)
      alone gives `phi >= 0`.  The best singleton is `C` in 29,019, `B` in
      17,897 and `A` in 9,542.
    - *Four factors* (`fm39/merge4.py`).  Merging both pairs gives
      `phi = sum_(d, d') phi(hat S_d hat S_d') - corr_pairing`, and each
      `hat S hat S` value is `T + W >= 0` by (E).  In 756 random words
      (`r <= 5`, labels `<= 9`), some pairing has `corr <= 0` in 491
      (65%), and `corr = 0` exactly in 72.  Assigned as FM-SEC19.
    - The other 8,792 words have every `M_Z > 0`, but small against
      `phi`, e.g. `M = 33` against `phi = 53,119` for `h_9 h_4^2 h_1^11`
      at `r = 2`.  The basis LP shows that a linear reduction to (E)
      forms is impossible for some of them.
  - **Theorem G0E4: four factors on their support region follow from (E)
    (FM-SEC19, luna_max_jupiter; verified by the main agent,
    `fm39/g0e4_check.py`).**
    - *Merge form.*  For a pairing `XY|ZV` at `e = 2r-4`,
      `phi = sum_(d in CG(X,Y), d' in CG(Z,V)) phi(hat S_d hat S_d' h_1^a) - corr`,
      with
      `corr = W(U_X U_Y U_Z, U_V) + W(U_X U_Y U_V, U_Z) + W(U_Z U_V U_X, U_Y)
      + W(U_Z U_V U_Y, U_X) - W(U_X U_Z, U_Y U_V) - W(U_X U_V, U_Y U_Z)`.
      Each `hat S hat S` value is `T + W >= 0` by (E).
    - *Support region.*  Let `lambda(I)` be the least fusion constituent.
      By the strict support rule `W = 0` for `p + q > N`, `corr` vanishes
      termwise exactly when:
      - `X + lambda(Y,Z,V) > N`, and the same for each singleton;
      - `|X-Z| + |Y-V| > N` and `|X-V| + |Y-Z| > N`.
      For sorted labels and the outer pairing `12|34` this is
      `s* = min_i (L_i + lambda(rest)) > N` and
      `L_1 + L_2 - L_3 - L_4 > N`.  That strictly contains Theorem LL4's
      region `mu_1 > N`, `mu_2 > N`; example `L = (6,5,4,3)`, `r = 2`,
      `a = 2`.
    - FM-CHK37 (luna_max_saturn, own exact code):
      - ACCEPT: the double merge (3,072 direct cases), LL4 containment
        with its strict gain, and the consumed (E) instances (`T + W` at
        pairs `(d, d')` with `d + d' <= N`).
      - REPAIR (applied): the strict support rule.  At equality a
        boundary `W` can be nonzero, e.g. `W(U_3,U_1) = 4` at `r = 2`,
        `a = 4`, or vanish by accident.
    - *Checks.*  1,366 words in the outer-pairing region (`r <= 5`,
      `a < 9`, labels `<= 11`): `corr = 0` exactly and `phi >= 0` in
      every one.  577 of them lie outside LL4.
    - The correction telescopes into phase-plane wedges, which gives an
      exact sign test but no universal sign rule.  Beyond four factors,
      iterated merges leave products of three or more `hat S`, which (E)
      alone does not control.  Those are the three-label plus words,
      FM-SEC23.
  - **FM-SEC24 (luna_max_pluto): candidates for a merge-closed
    induction (M5).**
    - Killed, with counterexamples:
      - insertion of `Q_t` with `|t| < 2` (`Q_0` at `r = 1` gives `-1`);
      - real-rootedness of the level generating polynomial
        (`1, 5, 35, 294` is not log-concave);
      - PSD Gram matrices of `hat S` insertions (`det = -1` for core
        `h_3 h_1`);
      - "some grouping has correction `<= 0`" (`h_2^3 hat S_3 h_1` at
        `r = 3`: every grouping has `C = 8`, `R = 38`).
      Each knob is strength that FM3 does not use.
    - Survivor, screened: the kernel-closed bound `|C_ij| <= R_ij` for
      every merge of two factors, with `Q_t` kernels (about 52,000 pair
      instances, no failure).
    - Main-agent classification: `R` and `C` are shared by the word with
      the two general factors and the word with the matching two
      `hat S` factors at the next lower level.  So `|C| <= R` holds
      exactly when FM3 holds for both of those words.  The survivor is a
      restatement of FM3 at the same label count, not a stronger
      induction hypothesis.  M5 still needs a hypothesis that controls
      `C` from data with fewer labels, such as a quantitative bound or a
      positive-form structure.
    - FM-SEC39 (luna_max_pluto): the 3-dimensional binary-form state is
      killed as an induction device (knob: an auxiliary state, not used by
      FM3).
      - Definiteness fails for 47 of 1,152 random consumer words,
        including the empty word.
      - Multiplication by a generator does not descend to `Sym_2`:
        `M_(0,2)(hat S_1 - hat S_3) = 0` while
        `M_(0,2)(hat S_1 (hat S_1 - hat S_3)) = (1,2,0)`.  So a
        merge-stable invariant needs a larger state whose kernel is
        invariant under multiplication.
      - The actual-vector half-space contained every sampled word, but
        that half-space is FM3 itself.
    - FM-SEC43 (luna_max_pluto): enlarged states.
      - Raw cross-kernel table `M_w(p,q)` with entrywise positivity:
        exact CG action, wrong cone (the empty word has entry `-1` at
        `(1,1)`, `r = 1`).
      - Finite level/suffix window: two words (`hat S_13`, `hat S_27`)
        with identical states are separated by one more generator.
      - GFM3 pointwise kernel cone: `Q_3 h_2 = -10` at `(0,0)`.
      - The E-polarized two-point state
        `lambda_w(P) = phi_(r + kappa(P))(w P)`, `P` in
        `{H_i H_j, H_i S_j, S_i S_j}`, carries every generator action
        exactly, but `H_2^3` has no finite positive expansion in pair
        probes (exact separator on labels `<= 6`, degree argument).  A
        closed unbounded completion is unsettled.
    - **FM-MECH28 (astra_max_ceres): candidate H, Hurwitz positivity of
      translation profiles.**
      - For a word `w` (`m` general `h`, any `hat S`, degree `d`), every
        `r >= max(1, ceil(m/2))` and `a >= 0`, put
        `P^w_(r,a)(t) = (1/2) E[(x-y)^(2r) (x+y)^a w(x+t, y+t)]`
        (suffix not translated).  Sign reversal gives
        `P = t^eps G(t^2)`, `eps = d + a mod 2`.
      - H: `G = 0`, or, after removing powers of `z`,
        `q(z) = b_0 z^n + ... + b_n` has `b_0 > 0` and every Hurwitz
        determinant `Delta_j(q) > 0`.  Since `Delta_n = b_n Delta_(n-1)`,
        H implies `phi_r(w h_1^a) >= 0`.  H holds for `w = 1`.
      - Exact insertion laws, uniform in the number of factors: with
        `(AP)_(r,a) = P_(r,a+1)`, `(RP)_(r,a) = P_(r+1,a)`,
        `L_t = A + 2tI` and
        `B(xi) = I - L_t xi + [2I + (L_t^2 - R)/4] xi^2 - L_t xi^3 + I xi^4`,
        `sum_k P^(w h_k) xi^k = B^(-1) P^w` and
        `sum_p P^(w hat S_p) xi^p = (2I - L_t xi + 2I xi^2) B^(-1) P^w`.
      - Screens: 3,904 profiles (single generators to label 40 and level
        32, words with up to 14 `hat S`, `r = 2` words `h h hat S hat S`),
        no failure; 420 insertion identities.  Rerun by the main agent.
      - Kills on the way (knob: strength FM3 does not use): endpoint
        deformations of `Q_t` (`phi_1(h_2 (4+xy)(4-xy)) = -1`), covariance
        (`-169884`), and every cone defined by nonnegative kernel moments
        (`F = (h_2 - 1)^2` is pointwise `>= 0` with nonnegative character
        coefficients, yet `phi_1(F hat S_4) = -1`; `F` is not a consumer
        word, and H excludes it).
      - Open: preservation.  The first step is `hat S_2` insertion,
        `Q = (1/2)[P_(r,a+2) + P_(r+1,a) - 4P_(r,a)] + 2t P_(r,a+1) +
        2t^2 P_(r,a)`.  H on single profiles cannot suffice (the array
        `P = 1` for even `a`, `t` for odd `a`, gives `4t^2 - 1`), so the
        proof must use compatibility of the actual arrays across `(r,a)`.
        FM-MECH29 (propagation) and FM-SEC49 (falsification) are on it.
    - **FM-MECH29 (astra_max_ceres): H is false** (knob: Hurwitz
      stability of the whole profile, while FM3 consumes only its constant
      term).
      - `w = hat S_2^9`, `r = 1`, `a = 0`:
        `P = 4q(t^2)`, where
        `q(z) = 128z^9 + 5184z^8 + 67968z^7 + 370944z^6 + 915264z^5 +
        1050768z^4 + 573720z^3 + 170928z^2 + 42858z + 8279`.
        The first seven Hurwitz determinants are positive, but
        `Delta_8(q) = -8938654984643632236257554566691734683648`.
        The consumer value is `phi_1(hat S_2^9) = 33116 > 0`.  The four
        parent profiles of `hat S_2^8` pass H, and the insertion law
        reproduces the child exactly.  Rerun by the main agent.
      - Survives: every actual profile array satisfies `C P = 0`, with
        `C = Pi(A,R) d/dt + 2T Pi(A,R) + 3 Pi_s(A,R)`,
        `Pi(s,D) = s^4 - 2(D+16)s^2 + (D-16)^2` (the boundary of the
        `(s, D)` support) and `(TP)_(r,a) = a P_(r,a-1)`.  It commutes with
        `A + 2tI` and `R`, so every insertion preserves it.  It excludes
        the artificial array `1 / t` but not `hat S_2^9`.
      - Killed on the way: PSD of joint `(s, D)` moment matrices (`det
        -16` for `hat S_2`, `r = 1`).
      - FM-SEC49 (falsification of H) was stopped as moot.
    - Large-level screen (main agent, `fm39/fm3_large_r_screen.py`).  At
      large `r` the weight concentrates at the antipodal corners
      `x = -y = +-2`, where every odd-label factor vanishes.  So words with
      several odd labels are decided by local forms there.  All 1,819
      words with up to 4 factors from `h_2..h_7`, `hat S_2..hat S_7`, at
      `r = 10, 20, 30, 40` and `a <= 3`: 29,104 exact evaluations, none
      negative.  (Per core, Theorem LR4 already leaves only finitely many
      `(r, a)`.)
    - FM-SEC82 (luna_max_uranus): level 1 as the exact coefficient
      inequality `b_00 + b_20 >= b_11` for the symmetric `SU(2)^2`
      expansion `F = sum b_mn U_m(x) U_n(y)` of the word.  The eigenvalue
      assignment sum needs cancellation, so a termwise route fails.  No
      proof.
    - **Randomized exact FM3 screen across levels (main agent,
      `fm39/fm3_random_screen.py`, fusion-kernel evaluator
      `fm39/fm3kern.py`).**  20,000 random consumer words with `r <= 8`,
      `m <= 2r` general `h` (labels `2..10`), `0..6` `hat S` (labels
      `2..8`), and suffix `a = 0..12`: 260,000 exact evaluations, none
      negative (130,738 zeros, least positive value 1).
    - **Unconstrained screen (main agent).**  The bound `m <= 2r` on the
      general `h` factors seems not to be needed.
      `phi_r(w h_1^a) = (1/2) E[(x-y)^(2r)(x+y)^a w]` is `>= 0` on all
      10,680 words with `r <= 3`, `2r < m <= 2r+3` general `h` (labels
      `2..4`), `0..2` `hat S` and `a <= 3`.  This matches the note's MP
      screens for general `kappa`.  A hypothesis without level
      bookkeeping (any `r >= 0`, any word) is therefore a live option for
      M5.  Genuineness alone is not enough: `F = (h_2 - 1)^2` is a
      genuine `SU(2) wr Z_2` character, yet `phi_1(F hat S_4) = -1`.
    - FM-SEC56 (luna_max_pluto), falsification of U-FM3.  No negative
      value on about 1,015,000 exact profiles: exhaustive `H_2..H_4`
      multisets with `m <= 2r+8` (`r <= 5`); up to four `H_2..H_8` with two
      `S_2..S_8`; pure `S`; `H_2^M` for `M <= 20`; random high-label words.
      Killed extensions (knobs not used by FM3):
      - real exponents: `E[|x-y| H_2] = -1024/(1575 pi^2) < 0` at
        `s = 1/2`, so the integrality of `2r` is essential;
      - pointwise nonnegative kernel factors:
        `E[H_2 Q_(-3)] = -1` at `r = 0`.
    - **Two exact frames for M5 (main agent).**
      - EVEN form, exactly FM3.  For labels `n_i >= 1` and signs `eps_i`
        with an even number of minus signs,
        `E[prod_i (U_(n_i)(x) + eps_i U_(n_i)(y))] >= 0`.  Minus labels 1
        are the level factors `x - y`, minus labels `L >= 3` are
        `h_(L-1)` times one `x - y`, and plus labels are `h_1` and
        `hat S_p`.  The value vanishes for an odd minus count, by the swap.
        It is multiaffine in `eps`, so FM3 is equivalent to positivity on
        the whole cube `[-1,1]^n`.  It is not real-stable (already
        `1 + eps_1 eps_2` at `n = 2`).
      - Super-dimension form (checked symbolically).  For `Spin(4)` with
        torus `(u_1, u_2)`: `x - y = ch S^+ - ch S^-` (the half-spin
        representations `V_1 (x) 1`, `1 (x) V_1`), and
        `(x-y)^2 = sum_k (-1)^k ch Lambda^k V`, where `V = V_1 (x) V_1`
        is the vector representation.  Hence
        `phi_r(w) = (1/2) sdim Inv_Spin(4)(M_w (x) Lambda(V (x) C^r))`
        for the genuine module `M_w` of the word.  In the EVEN form the
        value is `sdim Inv((x)_i W_i)`, with
        `W_i = V_(n_i) (x) 1 + Pi^[eps_i = -1] (1 (x) V_(n_i))`.  Candidate
        mechanisms: an odd differential with even cohomology, skew Howe
        duality on `Lambda(V (x) C^r)`, or a Lie superalgebra with even
        part `sl_2 + sl_2`.  FM-SEC58 (luna_max_mercury) is testing
        them.
      - FM-SEC58 (luna_max_mercury): exact skew Howe weights
        (`SO(4) x O(2r)`).
        - `s_r(a,b) = E[(x-y)^(2r) U_a(x) U_b(y)] = (-1)^b d_r(a,b)`, with
          `d_r(a,b) >= 0` the multiplicities in `(x+y)^(2r)`.  So
          `phi_r(w) = (1/2) E[(x+y)^(2r) chi_w(x,-y)]`.
        - The unbounded positive family it yields, products of `hat S`
          with even labels, is the known parity family.
        - Killed: a differential acting on the exterior factor alone (the
          channels have negative signed weight); a factorwise
          superalgebra extension with the prescribed grading.
        - The Howe formula is exact but not termwise positive
          (`h_2^4`, `r = 2`: `2894 - 2880`).  No movement on the full cone.
    - **FM-MECH30 (astra_max_ceres): Hypothesis B, the leading M5
      candidate** (reproducer `fm39/mech30_hypB_repro.py`, rerun
      exactly by the main agent).
      - For a label list `lambda = (n_1, ..., n_L)` put
        `m(S) = [U_0] prod_(i in S) U_(n_i)`, and on
        `G_L = F_2^L / <(1,...,1)>` put `f_lambda([S]) = m(S) m(S^c)`.
      - B: `f_lambda` is a nonnegative combination of indicators of
        linear subspaces of `G_L`.
      - B implies the whole EVEN statement, i.e. FM3 for every consumer
        word at every level with any number of factors.  At an even
        minus set `T`, the value is `2 sum_(W <= ker chi_T) a_W |W| >= 0`.
      - Equivalently (main agent): the sign-group action on the invariant
        space `Inv (x)_i (V_(n_i) (x) 1 + 1 (x) V_(n_i))` is a positive
        rational combination of permutation representations.  So a basis
        of invariants permuted by the sign flips would prove FM3.
      - Screens: exact certificates (LP, then exact rational
        reconstruction) for all sorted lists of length 2..7 with labels
        `<= 4` (164), length 8 with labels `<= 3` (25), and
        `(1,1,2^k)` for `k <= 12` (via binary-code weight enumerators).
        The last family includes `hat S_2^9`, which killed H.
      - Proved: merger restriction,
        `f_(mu,a,b)|_(bit_a = bit_b) = sum_(d in CG(a,b)) f_(mu,d)`.
      - Delimiting kills: (i) abstract gluing fails.  On `F_2^4`, a
        function with nonnegative Fourier transform whose every
        hyperplane restriction is in the subspace cone is not itself in
        the cone.  So a proof must use the actual fusion tables.  (ii)
        Intermediate-spin cutoffs `F_J` break Fourier positivity
        (`(1,1,2)`, `J = 1`: value `-1`).
      - Open: the insertion step `(Ext)`.  FM-MECH31 seeks a canonical
        decomposition or a permuted basis.  FM-SEC68 (luna_max_eris)
        tries to break B at scale.
      - Necessary condition, main agent (`fm39/hypB_necessary.py`).
        Every subgroup indicator satisfies
        `1_K(T_1+T_2) + 1_K(0) >= 1_K(T_1) + 1_K(T_2)`.  So B implies
        `F(T_1+T_2) + F(0) >= F(T_1) + F(T_2)` for the sign-pattern values
        `F`.  It holds on 12 random lists of length 6..9 with labels
        `<= 7`, beyond FM-MECH30's range.
      - All-ones lists: pairs of noncrossing matchings on `S` and `S^c`
        give subspace indicators when their union is noncrossing, but
        affine cosets (2-colourings of the crossing graph) when it is not.
        So a canonical decomposition must re-pair crossing configurations.
      - All-ones, `L = 8`, explicit (main agent):
        `f = 4 * 1_E + (1/6) sum_(|T| = 3) 1_(W_T) + (2/3) 1_0`, where
        `W_T` is the span of the pairs inside the triple `T`.
      - **FM-SEC68 (luna_max_eris): Hypothesis B is FALSE** (knob:
        subspace-mixture membership, stronger than the Fourier positivity
        FM3 uses; `fm39/sec68_hypB_counterexample.py`, rerun exactly).
        - `lambda = (1^8, 6)`.  With representatives having the label-6
          bit zero, `f = 7 delta_0 + sum_(|x| = 2) delta_x` on `F_2^8`.
        - Main-agent proof.  A subspace in the support of `f` has only
          weight-2 nonzero vectors, so it is `{0, v}` or a triangle
          `{0, e_ij, e_jk, e_ik}`.  Every subspace contributes to the
          origin, whose budget is 7, and covers at most 3 of the 28
          weight-2 vectors.  A decomposition would need total
          coefficient `>= 28/3 > 7`.
        - Eris's separator: `y` = `(3,-1,-1,3,3,-1,-1,3,3)` by weight has
          `sum_(W) y >= 0` on all 417,199 subspaces, but
          `sum y f = -7`.
        - FM3 holds at this list: `F(u) = 3 + 2(4 - |u|)^2 >= 3`.
        - So the sign-group action on the invariant space is not a
          positive combination of permutation representations in general.
          FM-MECH31's theorem for repeated odd labels stands, and FM-CHK42
          accepted it.
      - **FM-MECH33 (astra_max_ceres): hypothesis H_AC, the new leading
        M5 candidate** (`fm39/mech33_hAC_repro.py`, rerun exactly).
        - No-go: any construction from subspace indicators by
          nonnegative sums, tensor and pointwise products, linear pullbacks
          and linearly constrained hidden variables stays inside B.  The
          family `(1^N, N-2)`, with `f = (N-1) delta_0 + sum_(|S|=2)
          delta_S`, kills all of them for `N >= 7`, including block-tensor
          relaxations.
        - H_AC: `f_lambda = sum_nu a_nu (p_nu * p_nu)` with `a_nu >= 0`,
          `p_nu >= 0`, and `(p*q)(S) = sum_X p(X) q(S+X)`.  Equivalently,
          the matrix `C(X,Y) = f(X+Y)` is completely positive,
          `C = V V^T` with `V >= 0` entrywise.  Plain PSD of `C` is FM3
          itself.
        - H_AC implies FM3, since the Fourier transform of `p*p` is
          `|p_hat|^2`.  Normalization (FM-CHK45 repair): with the
          transform over `G = F_2^L/<1>`, the EVEN value is
          `2 f_hat_G(T)`.  It is strictly between Fourier positivity and B:
          `1_W = |W|^(-1) 1_W * 1_W`, and the gluing control `g` on
          `F_2^4` fails it.
        - At `(1^8, 6)`: `f = (1/2) p*p + 3 delta_0`, with `p` the weight-1
          indicator.
        - Proved families:
          - inserting the label `n = sum(mu) - 2` into any list `mu`:
            `f = (1/2) p*p + (L - 1 - t/2) delta_0`;
          - `(1^N, n)` with `3n >= N-2`, by sphere autocorrelations and an
            explicit binomial identity.
        - **FM-MECH34 (astra_max_ceres; `fm39/mech34_hAC_d2_repro.py`,
          rerun exactly): H_AC for every list with a label
          `n = sum(mu) - 4`**, for any number of mixed labels.
          - The table is
            `f = [C(L,2) - t] delta_0 + (L-3) E_11 + E_22 + E_112 +
            2 E_1111`, where `t, v` count the labels 1 and 2 in `mu` and
            `E_rho` indicates the subsets with label pattern `rho`.
          - The factors are `E` (single 1's), `A` (pairs of 1's), `B`
            (single 2's) and `A` plus weighted 2-points.  The certificate
            is explicit in `(t, v, w)`, with every coefficient `>= 0`.
          - Checked on 87,092 direct fusion entries in 400 lists.  With
            FM-MECH33 this gives H_AC, hence FM3 at every level, for every
            list with `2 max(lambda) >= sum(lambda) - 4`.
          - Its factors are pointwise disjoint in the inserted bit, so they
            also give the compatible insertion data for
            `n = sum(mu) + h - 4`.
          - Boundary screens: `(1^10,2)` and `(1^12,2)`; all 21 inner-range
            `(1^N, n)` with `N <= 20` and `3n < N-2` (two-sphere atoms);
            `(1^(2m),2,2)` for `m <= 6` (cycle-code atoms for `m = 5, 6`).
          - Next distance: at `n = sum(mu) - 6` the table (checked on
            22,016 entries) has six-strand patterns `E_111111`, `E_11112`,
            `E_1122`, `E_222`, `E_1113`, `E_123`, `E_33`, and no factor
            rule is known yet.  Since every list has `sum(mu) - n = 2d`
            with `n` its largest label, a construction uniform in `d` is
            H_AC for the full cone.
      - FM-SEC76 (luna_max_eris): no FM3 counterexample in large-label
        screens.
        - Two hats: 53,760 full-suffix cases (labels `<= 20`, `r <= 4`) and
          36,462 boundary-band cases (labels `<= 60`).
        - Three hats `hat S_2^2 hat S_p` (`p <= 60`): 25,080 cases.  Also
          `h`-cores with one or two factors to label 40.
        - Tight normalized families (value 1 at `r = 1`):
          `hat S_1 hat S_q` and `hat S_2 hat S_q` (rate `q^4 4^(-q)`),
          `hat S_2^2 hat S_p`, `h_n`, `h_2 h_n`, `h_3 h_n`.  For
          `hat S_n^2` the minimum is `3/(4(n+1)^2)`, at `a = 0`.
        - Its first unresolved family,
          `phi_1(hat S_p hat S_(p+1) h_1^(2p+1))`, is already proved
          (FM-SEC78, accepted by FM-CHK44).
      - FM-SEC88 (luna_max_venus; `fm39/sec88_hAC_12sector_repro.py`,
        rerun exactly): H_AC on the `{1,2}` sector.
        - The profile is `m(i,j) = Delta^j Cat_(i/2)` for even `i`.
        - Bi-radial certificates exist for all even `N` with
          `N + M <= 6` and for `(0,7), (2,5), (6,1)`.  They are exactly
          ruled out at `(1^4, 2^3)` (separator `y . f = -1`), but full
          H_AC holds there with a 29-subgroup (B-type) certificate.
        - Subgroup-square certificates for `(1^10, 2)` and `(1^12, 2)`,
          and a bi-radial certificate for `(2^10)`.
        - No uniform rule for the sector.
      - Literature check (main agent).  FM3 is Ginibre's single-site
        condition Q3 for the cone of SU(2) characters (the O(4) zonal
        functions `U_n`).
        - Sylvester, "The Ginibre inequality", Commun. Math. Phys. 73
          (1980) 105-114, shows that the multi-site Ginibre inequality
          (positive definiteness of certain functions on the cycle group
          of a graph) fails on some graphs for spin dimension `>= 3`.
        - Herbst, arXiv:2209.11850, proves the Griffiths inequalities
          for non-interacting rotors; Tokushige, J. Stat. Phys. 192
          (2025) 25, gives a graphical proof for the XY model.
        - The known failures concern O(N) interactions `sigma_i . sigma_j`
          (matrix coefficients) on graphs with cycles.  The class-function
          Q3 here does not feed into them.  No conflict with FM3 is known,
          and no result on the SU(2) character cone was found.
      - FM-CHK48 (luna_max_neptune, own code) on FM-MECH35: ACCEPT all
        five items.
        - The `d = 3` table is rederived and checked on 520 fusion
          profiles.
        - The `t <= 3` residual formulas and bounds are checked (exact
          residual scan to `L = 100`).
        - The `v <= 1` identity and budget polynomials are checked (279
          profiles), as are the small cases.
        - The checker's item-5 REPAIR (a `(q^5)` formula "false at
          `q = 2`") is withdrawn.  It evaluated the six-label list
          `(2^5, 4)`, where the table is `(15, 2, 1)`.  For `lambda = (2^5)`
          the main agent confirms the stated `(6, 1, 1)`.
      - FM-CHK47 (luna_max_eris, own code, no imports of the
        reproducers).
        - ACCEPT: Claim R, including the Yang-Baxter check, the norm, the
          one-colour q-Hermite recursion and the Wick positivity.
        - ACCEPT: Claim Z, on 240 cases.
        - ACCEPT: the component-model identity, on 1,088 subset counts.
        - ACCEPT: the Fock evaluator (280 cases) and both Bernstein
          screens.
        - ACCEPT: the graph-product screens at `N = 2, 3`.
        - REPAIR: the product-state convention in the cumulant reading,
          and the sampled global minimum of the crossing-set kill (both
          applied).  Also the self-test import in `qs_fock.py` (fixed).
      - FM-CHK46 (luna_max_saturn, own code): ACCEPT FM-MECH34 Theorem 1
        (21,975 fusion cases; the `rho` minima `17/12, 5/3, 17/12`;
        adjacency-free supports), the FM-SEC77 `q = 1` rotation theorem
        (Hermite addition via the generating function; 300 random lists),
        and the FM-SEC79 certificates (`a = 5, 6` against direct values).
      - Level 1 as an operator statement (main agent).  For
        `R = prod_k (V_(n_k) (x) 1 + 1 (x) V_(n_k))`, level 1 for every
        word is `R[a][b] <= sum_(c in a (x) b) R[c][0]` for all
        `a, b >= 1`.  Equivalently, `M_rho - K_R` has nonnegative entries
        in the character basis, with `K_R` the integral operator with
        kernel `R` and `rho = int R dh`.  Multiplying `R` by `S_n` gives
        `M_(chi_n)(M_rho - K_R)`, which is `>= 0`, plus a term with
        entries `sum_(c in a(x)b) R[c][n] - sum_(d in n(x)b) R[d][a]`.
        These are antisymmetric under `a <-> n`, so a termwise induction
        fails; this is the cross-term obstruction again.
      - FM-CHK45 (luna_max_mercury, own code) on FM-MECH33:
        - ACCEPT: `B` is contained in H_AC, and the gluing control is
          outside it (support argument, origin mass `3 < 5`).
        - ACCEPT: the B-class no-go, with `C(N,2) <= 3(N-1)` failing for
          `N >= 7`.
        - ACCEPT: the boundary insertion `(L-1) delta_0 + sum delta_(e_i+e_j)
          = (1/2) p*p + (L-1-t/2) delta_0` (300 lists).
        - ACCEPT: the sphere certificate `f_(1^N,n) = sum_h b_h
          (p_(k-h) * p_(k-h))`.  Only `b_1 = (n-k+1)/(k(k+1))` can be
          negative, and it is `>= 0` exactly when `3n >= N-2`.
        - REPAIR: the Fourier normalization (applied above).
      - FM-SEC89 (luna_max_mars; `fm39/sec89_hAC_census_repro.py`,
        assert-only, reruns exactly): structured dictionaries.
        - Class-profile Hamming spheres fail at `(1,1,1,1,2,2)`, mask 19.
          Fusion-path factors `P_(j,k)(S) = m_j(S) m_k(S^c) +
          m_k(S) m_j(S^c)` repair it, with coefficients `2, 1/2, 1/4` on
          `P_(0,8), P_(1,7), P_(2,6)`.
        - At `(2^6, 4)` spheres, path factors and point pairs all fail,
          with an integer Farkas vector.  H_AC still holds there:
          `f = 15 delta_0 + (1/8) sum_(O_1) 1_H*1_H +
          (1/24) sum_(O_2) 1_H*1_H`, with subspace orbits of sizes 15 and
          30.
        - The gap-two insertion law is rechecked on 200 lists.
        - No H_AC counterexample.  The CP case `(1^13, 3)` flagged by
          FM-SEC85 is among the 21 inner-range lists that FM-MECH34
          certified with two-sphere atoms.
      - **FM-MECH35 (astra_max_ceres): H_AC at distance `d = 3` for
        `t <= 3` or `v <= 1`** (`t, v` the counts of labels 1 and 2 in
        `mu`).  Reproducers were extracted from the run log and rerun
        exactly: `fm39/mech35_hAC_d3_repro.py` (1,012 lists, 519,830
        entries), `fm39/mech35_hAC_d3_budget.py` and
        `fm39/mech35_hAC_d3_witnesses.py`.
        - The `d = 3` table is
          `F_6 = D delta_0 + [C(L-2,2)-t+2] E_11 + (L-3) E_22 +
          (L-4) E_112 + 2(L-5) E_1111 + E_33 + E_123 + E_1113 + E_222 +
          2 E_1122 + 3 E_11112 + 5 E_111111`, with
          `D = C(L+1,3) - t(L-1) - v`.
        - For `t <= 3` the factors are built from `E, A, H` (singles,
          pairs, triples of 1's), `B` (single 2's), `K` (single 3's),
          `J = E*B`, and `R = (1/2) sum Q(delta_i + delta_(j+k))` over
          triples of 2's.
        - For `v <= 1` the key identity is
          `Q(5H+6J+12K_0)/144 + Q(H+3J)/18 + Q(H+4K_0)/48`.
        - The origin budgets are proved by polynomials with positive
          coefficients.  The small cases (36 count patterns) are certified
          exactly.
        - Two restricted rules are killed, and H_AC holds at both
          witnesses:
          - factors supported on whole-cluster subsets of weight
            `<= d + c`: `(2^4; n = 2)`, and `(q^5)` needs support radius
            `>= 4d/3`;
          - uniform weight-class factors: the Gram determinant is `-1/64`
            at `(1^6, 2^2; n = 4)`.
        - Open at `d = 3`: `(1^t, 2^v)` with `t >= 4`, `v >= 2`.  No rule
          uniform in `d` yet.  Per `d` this is becoming a sequence of
          separate certificates, so the ladder is paused in favour of an
          insertion-stable hypothesis (FM-MECH36).
      - FM-SEC113 (luna_max_venus; `fm39/sec113_hACq_L8_repro.py`, rerun
        exactly):
        extended H_AC_q screen.
        - All 25 even-total lists of length 8 with labels `<= 3`
          (including `(1,1,2,2,3,3,3,3)`, unresolved in FM-SEC101) have
          coefficientwise certificates.  The dictionary is subgroup-orbit
          squares, q-independent: 775 coefficient vectors, 2,844 terms.
        - Structured lists of length 9-11 are certified:
          `(1^8,2), (1^6,2^3), (1^10), (1^9,3), (1^10,4)`.
        - The q-graded Fourier transforms are positive on 1,250,432
          entries (lengths 8-10); the minimum is 2.
        - `(1^8, 2^2)` at `q^4` lies outside the tested finite dictionary
          (exact dictionary separator).  Its Fourier minimum is 814, so
          it is not a counterexample.
      - **FM-MECH38 (astra_max_ceres; `fm39/mech38_hACq_d2_repro.py`,
        rerun exactly): H_AC_q PROVED for `d <= 2`, uniformly**, i.e.
        for every list with `2 max(lambda) >= sum(lambda) - 4`.
        - Lemma: the coefficients of `H_(N-2)` and `H_(N-4)` in
          `prod H_(mu_i)` are the q-positive sums
          `J(mu) = sum_i [s_i][mu_i]` and
          `K(mu) = sum_i [s_i choose 2]_q [mu_i choose 2]_q [2]! +
          sum_(i<j) [s_i][mu_i][s_j - 2][mu_j]`, with `s_i` the prefix
          sums.
        - `f_q = [n]! g_q`.  At `d = 1`, `g_q = (1/2) E*E + (J - t/2)
          delta_0`.  At `d = 2`, `g_q = G_0 + q G_1 +
          (1/2) J_(>=2)(q) E*E + R_(>=2)(q) delta_0`, with explicit
          q-independent factors and `R = K - (t/2)J` coefficientwise
          nonnegative, proved by weighted q-integer comparisons.
        - KILL: factorizing each aggregate fusion trajectory separately.
          It fails at `(1^4)` with Fourier `q - 1`.  After averaging
          over block orders it still fails, at `(1,1,1,1,2)` with
          `-1/15`; 343 of 45,961 averaged tables fail.
        - A minimal cross-trajectory transfer `delta_0/15` repairs that
          witness exactly.
        - Knob: trajectory separation.  A general rule must combine
          trajectories within a q-degree.
        - FM-CHK51 (luna_max_eris, fresh code): ACCEPT all four items
          (3,095 ordered lists with labels `<= 4`, total `<= 12`; `d = 1`
          on 496,656 entries, `d = 2` on 496,584).
      - FM-SEC119 (luna_max_mercury; `fm39/sec119_d3_radial_repro.py`,
        rerun exactly): distance three, two or more labels 1.
        - The `d = 3` table is `g_q(S) = m_q(mu_S) K_(3-k)(mu_(S^c))`.  The
          low-side moments are listed by label type, and the append
          recurrence is
          `K_d(nu a) = sum_j K_j(nu) [S-2j choose d-j]_q
          [a choose d-j]_q [d-j]_q!`.
        - KILL: radial factors `aE + bH`, `A` on the 1-labels.  They force
          `g(4) >= (9/10) g(6)`, which fails at `(1^7)`, `n = 1`
          (`4 < 9/2`).  This kills only the radial dictionary: that list
          is B via non-radial subgroup orbits (FM-MECH31, FM-SEC113).
        - Comparison table of the `d = 1, 2, 3` factor families recorded.
      - Main-agent screen (`fm39/twoM_block_codes_screen.py`): for
        `(2^M)`, a radial B certificate over block-indicator codes plus
        even-weight codes exists only for `M <= 6`, and is infeasible for
        `M = 7..14`.  General subgroup orbits certify `M <= 10`
        (FM-SEC88, FM-SEC113).  So the needed codes grow in complexity
        with `M`.
      - **FM-MECH41 (astra_max_ceres; `fm39/mech41_B_12lists_repro.py`, all
        assertions pass): B holds for every list with labels in {1, 2}.
        So `phi_r(h_1^a hat S_2^k) >= 0` for all `r >= 1`, `a, k >= 0`.**
        - The {1,2} sector is closed at every level with any number of
          factors, including the odd-`k`, even-`a` branch left open by
          FM-SEC69/72 and the repeated-2 sector.
        - Lemma: the semicircle quantile `v(t)` has
          `v(t) = sum_k A_k (pi t/2)^(2k+1)` with `A_k > 0`, by Lagrange
          inversion of `int_0^w sqrt(1 - s^2/4) ds = w(1 - rho(w^2))`;
          convergence on `[0,1]` by Picard iteration.
        - Positive Walsh realization.  Take independent signs
          `eps, sigma_1, sigma_2, ...`, `u = 1/2 + sum 2^(-j-1) sigma_j`
          (uniform) and `X = eps v(u)` (semicircle).  Then `X` and
          `U_2(X) = X^2 - 1` have nonnegative Walsh-Fourier coefficients
          `a_1, a_2` on the Boolean group; `E[X^2] = 1` removes the
          constant.
        - Hence `f(S) = m(S) m(S^c) = sum over charge assignments gamma`
          with `xor gamma_i = 0` of `(prod_i a_(n_i)(gamma_i))
          1_(W_gamma)(S)`, where `W_gamma = {S : xor_(i in S) gamma_i = 0}`
          is a subgroup.  This is a nonnegative subgroup mixture, hence
          B, H_AC and FM3.  The construction is global (no insertion),
          so FM-MECH40's obstruction does not apply.
        - The canonical extension to other labels fails: `U_3(X)` has
          Walsh coefficient `-(2/5) E|X|` at the sign `eps`, and the grouped
          `E_4` coefficient at `(1,1,1,3)` is `< -1/25`.  That table is
          still B, by a three-plane regrouping.
        - FM-CHK53 (luna_max_eris, fresh code): ACCEPT all four items.
          - Lemma 1: `A_k > 0` for `k <= 30`, and the inversion is checked
            through degree 61.
          - The Walsh construction and charge identity, including a
            moment-matching finite model and 35 words with 1,792 entries.
          - The consumer translation.
          - The `(1,1,1,3)` bound and the three-plane repair.
        - Main-agent observation: B is known to fail only on lists with a
          label `>= 5` (`(1^N, N-2)` for `N >= 7`; `(1^8,6)`), and the
          `L = 8` certificates for labels `<= 3` (FM-SEC113) are all
          subgroup sums.  A realization with `U_1..U_4` Walsh-nonnegative
          would prove FM3 for every word with labels `<= 4`.  This is
          FM-MECH42.
        - Lemma (main agent; checked numerically on 300 random
          coefficient vectors).  Let `X = sum a(gamma) chi_gamma` with
          `a >= 0` and `sum a^2 = E X^2 = 1`.
          - `U_2(X) = X^2 - 1` is automatically Walsh-nonnegative.
          - CORRECTION (FM-MECH42): `U_4(X)` is not automatic.  The
            earlier 'three pair-splittings without overlap' count was
            wrong.  There are six families (an equal pair plus a
            complementary pair summing to `gamma`), and they overlap on
            triples, giving only
            `U_4 coefficient >= b(gamma)(3 - 8 max a^2)`.  The FM-MECH41
            realization has a negative `U_4` coefficient at `sigma_1`.  The
            numerical check had missed this (random vectors with small
            max coefficient).
          - `U_3(X)` has coefficients `a^(*3) - 2a >= a(3 - 2a^2) - 2a`,
            which is nonnegative wherever `a <= 1/sqrt 2`.
          - The quantile realization fails at `U_3` only because
            `sgn(X)` is a single character (the top dyadic digit), with
            coefficient `E|X| = 8/(3 pi) > 1/sqrt 2`.  Every dyadic encoding
            has this defect.
          - (Superseded by FM-MECH42: no such realization exists.)
          - `U_5` and `U_6` cannot be automatic: B fails at `(1^7,5)` and
            `(1^8,6)`.
          - Necessary condition, depending only on the law: if `sgn(X)`
            equals a character `chi`, the `U_3` coefficient at `chi` is
            `E|X|^3 - 2 E|X| = -16/(15 pi) < 0` for every realization.  So
            a `U_3`-nonnegative realization needs `sgn(X)` not to be a
            character.
          - For three signs, all 168 nonnegative-Walsh bijections onto the
            dyadic grid have a character as sign.  So 'uniform, then
            quantile' looks closed.
          - Candidate shape: `X = G(Z)` with `Z = sum_j c_j sigma_j`.
            `sgn(Z)` is majority-like and each single-sign coefficient
            `c_j` is small.  `G` must be odd with positive Taylor
            coefficients on the range of `Z`, hence convex on
            `[0, max Z]`.
      - **FM-MECH42 (astra_max_ceres;
        `fm39/mech42_realization_nogo_repro.py`, rerun): the
        positive-Walsh realization route stops at {1,2}.**
        - Theorem: every semicircle realization on a Boolean group with
          nonnegative Walsh coefficients `a` has exactly one character
          `gamma*` with `a(gamma*)^2 = p > 1/2`.  Its `U_3` coefficient
          satisfies `b(gamma*) <= -(1-p)(2p-1)/sqrt p < 0`, and all other
          `U_3` coefficients are `>= 0`.
        - Proof: `b = a(1 - 2a^2) + R` with `R >= 0`.  Parseval gives
          `sum a b = E[X U_3(X)] = 0`, hence `sum a^4 >= 1/2`.  If equality
          held, `X` would be `(chi + eta)/sqrt 2`, which has `E X^6 = 4`, not
          5.  Main-agent check of the proof: correct.
        - So no common realization makes `U_1` and `U_3` both
          Walsh-positive, whatever the encoding.  `U_4` is not automatic
          either: `E[U_4(X) sigma_1] < 0` for the FM-MECH41 realization.
        - `(1,1,1,3)` is still B (`f = 1_<1111>`).
        - Knob: one common realization for all labels.  B and H_AC allow
          list-dependent decompositions, which FM3 needs anyway.
        - Necessary condition (main agent).  In any realization where
          `U_n` and `U_m` (`n != m`) both have nonnegative Walsh
          coefficients, their supports are disjoint, because
          `sum_gamma a_n(gamma) a_m(gamma) = E[U_n U_m] = 0`.
          - For {1,2} this holds through the sign-charge parity.
          - For {2,4}, with `Y = U_2(X)` having coefficients `c`, it forces
            `c*c = c` on `supp(c)`.  That is consistent with `E Y^3 = 1`
            and `E Y^4 = 3`, but excludes a single subgroup
            (`|H| = (5 +- sqrt 5)/2`).  FM-SEC123 is mapping the
            admissible label sets.
      - **FM-MECH44 (astra_max_ceres; `fm39/mech44_label3_repro.py`, rerun
        exactly): the {1,2} sector plus one label 3 is PROVED.
        `phi_r(h_2 h_1^a hat S_2^b) >= 0` and
        `phi_r(hat S_3 h_1^a hat S_2^b) >= 0` for every `r >= 1`,
        `a, b >= 0`.**
        - Theorem: `M_(b+1) >= |Q_b|` for all even `A, E`.
        - Radial coordinates: `d = t s`, `z = s^2 (1+t)^2/16`, so that
          `Z = beta z - 2` and `Z - P = alpha z - 2`, with
          `beta = 8(1+t^2)/(1+t)^2 >= 4` and
          `alpha = 4(1+3t^2)/(1+t)^2 >= 3`.
        - Then `M_(b+1) - Q_b = (16^(N+1)/pi^2) int t^E/(1+t)^(2N+2)
          int_0^1 z^N sqrt((1-z)(1-kappa z)) (beta z - 2)^b (alpha z - 2)
          dz dt`.
        - Every inner integral is `>= 0` by a positive radial IBP
          identity (Lemma 3) together with a reflected-weight lemma for
          `(beta z - 2)^k` (Lemma 2).  The small cases use genuine
          characters.
        - General extra label: the IBP recurrence for
          `R_(k,b) = E[s^A d^E P^k Z^b]` eliminates every
          `phi_r(W h_1^a hat S_2^b)`, `W in {h_(n-1), hat S_n}`, into an
          explicit combination `sum_j c_j M_j` (uniform in `n`; 1,200
          direct checks).  Its positivity is open.
        - KILL: raywise positivity for the general label.  `h_3 h_1^3` has
          ray value `-1/11` at `t = 1` (true value 1), and
          `hat S_6 hat S_2` at `r = 2` has `-36/35` (true value 1).
        - Next requirements: for `h_3`, `M_(b+1) - M_b >= 0` (A odd-shifted);
          for `hat S_4`, `M_(b+2) + M_(b+1) - 2 R_(2,b) >= 0`.
        - FM-CHK54 (luna_max_eris, fresh code): ACCEPT all four items.
          - The identities and the consumer translation.
          - Lemmas 2 and 3.
          - The radial coordinates, Jacobian and normalization, and the
            small cases (576 parameter checks).
          - The recurrence (1,470 cases), the elimination (1,512), the
            closed forms for `n = 3..12`, and both ray counterexamples.
      - **FM-MECH45 (astra_max_ceres; `fm39/mech45_extra_label_456_repro.py`,
        rerun exactly): the {1,2} sector plus one label 4, 5 or 6 is
        PROVED.  `phi_r(h_(n-1) h_1^a hat S_2^b) >= 0` and
        `phi_r(hat S_n h_1^a hat S_2^b) >= 0` for `n = 4, 5, 6`, every
        `r >= 1`, `a, b >= 0`.**
        - Coordinates: `x = 2 sqrt z`, `y = 2 q sqrt z` (fold to
          `x >= |y|`), `f = beta z - 2 = Z`, `P = c (f+2)`,
          `c = q/(1+q^2)`, `beta = 4(1+q^2)`, `B = beta - 2`.  Write
          `W = s^eps F(Z,P)`, `A = a + eps`, `E = 2r`, `N = (A+E)/2`.  The
          consumer is `4^(N+1)/pi^2 int (1+q)^A (1-q)^E I_b(F) dq` with
          `I_b(F) = int_0^1 z^N sqrt((1-z)(1-q^2 z)) f^b F dz`.
        - Lemma 2 (radial moment bounds, `N >= 2`):
          `0 <= B^j J_b - J_(b+j) <= 2 beta j B^(j-1) J_b/(N+b+3)`.  It
          follows from the reflected-weight lemma and the radial IBP
          identity of FM-MECH44.
        - Proposition 3 (criterion for any polynomial `F`): if
          `(N+b+3) F(q,B) >= 2 beta D(q)` on `[-1,1]`, with
          `D >= sum_j j max(a_j,0) B^(j-1)`, then `I_b(F) >= 0` on that ray.
        - Cutoffs `K` (criterion holds for `N + b >= K`, by exact degree-64
          Bernstein certificates): `h_3` 6, `hat S_4` 12, `h_4` and
          `hat S_5` 24, `h_5` 25, `hat S_6` 51.  `N = 1` is positive for
          every label (a character argument).  The finite remainder
          `2 <= N < K`, `b < K - N` is 28,946 exact Catalan evaluations:
          all `>= 0`, with 10 zeros.
        - Proposition 7: for every label `n` a finite cutoff `K_n^+-`
          exists, because `F_n^+-(q,B) > 0` on `[-1,1]`.  So each fixed
          label is proved outside a finite box.  The box grows with `n`
          (`K - n/2` = 4, 10, 21.5, 22, 48 above), so this is not uniform
          in `n`.
        - Proposition 8 (KILL, raywise positivity inside the box): at
          label 7, `phi_1(h_6 h_1^2 hat S_2^3)` and
          `phi_1(hat S_7 h_1 hat S_2^3)` have ray values `-32/105` at
          `q = -+1` and true value 4.  Knob: positivity before angular
          integration.  Main-agent remark: both lists are at distance
          `d = 1`, inside the proved H_AC_q region (FM-MECH38).
        - Main-agent reformulation, uniform in the label.  By the `x <-> y`
          symmetry, `phi_r(hat S_n h_1^a hat S_2^b) = [U_n] g_(2r,a,b)` and
          `phi_r(h_(n-1) h_1^a hat S_2^b) = [U_n] g_(2r-1,a,b)` (384 direct
          agreements, `n <= 8`), where
          `g_(k,a,b)(x) = E_y[(x-y)^k (x+y)^a (x^2+y^2-2)^b]`.  So the whole
          one-extra-label sector, for all `n` at once, is:
          every `g_(k,a,b)` is a nonnegative combination of the `U_n`.
          `fm39/extra_label_upositivity_screen.py`: 822 cases with
          `k + a + 2b <= 24`, no negative coefficient.  Only 68 of them are
          monomial-positive, so the ballot expansion of `x^m` alone does
          not suffice.
        - Also: the value vanishes unless `n <= 2(N+b) - eps`.  The
          distance of the list is `d = N + b - (n+eps)/2` (`hat S_n`) or
          `N + b - (n+1+eps)/2` (`h_(n-1)`), so only `d >= 3` is open.
        - Main-agent exact formula (`fm39/extra_label_determinant_formula.py`,
          9,790 entries with `k + a + 2b <= 22`, no mismatch).  With half
          angles `X = p rho`, `Y = p/rho`:
          `x - y = (rho - 1/rho)(p - 1/p)`, `x + y = (rho + 1/rho)(p + 1/p)`
          and `Z = (p^2 + p^-2)(rho^2 + rho^-2) + 2`.  So
          `d^k s^a Z^b = sum_i C(b,i) 2^(b-i) A_i(p) A_i(rho)` with
          `A_i(t) = (t - 1/t)^k (t + 1/t)^a (t^2 + t^-2)^i`, and the Weyl
          constant-term formula gives
          `[U_n] g_(k,a,b) = sum_i C(b,i) 2^(b-i) (u_n v_(n+2) - u_(n+2) v_n)`,
          where `u_j = [t^j] A_i` and `v_j = u_(j-2) + u_(j+2)`.
          - Single-`i` terms can be negative (2,763 of 9,790), so for `b > 0`
            the sum over `i` is needed.
          - `b = 0`: with `alpha_m = [T^m](T-1)^k (T+1)^a` (a Krawtchouk
            value), `D_m = alpha_m^2 - alpha_(m-1) alpha_(m+1)` and
            `M = (k+a+n)/2`, the formula reads `[U_n] g_(k,a,0) = D_M - D_(M+1)`.
            So `b = 0` says exactly that these Turan determinants decrease
            away from the centre `m = L/2`.  They are `>= 0` by Newton's
            inequalities.  The three-term recurrence alone does not give
            the decrease: the quadratic form in `(alpha_(m-1), alpha_m)` is
            indefinite.
          - CORRECTION (main agent): both forms are already in the note.
            `b = 0` is Lemma FM27, proved at every level by Theorem OL
            (FM-MECH6).  The `i`-rows are the cyclotomic-multiple rows
            `(1 + z^2)^i P` of the (E) analysis.  At `b = 1` the uniform-`n`
            statement is exactly the two-label inequality (E) at `q = 2`:
            the minus sign is proved (OL on `(1 - z^2) P`), and the plus sign
            is the open FM-SEC35/38/44 item.  So the {1,2} sector plus one
            arbitrary label contains that open item.
          - Corollary of FM-MECH41/44/45 (main agent).  (E) at `q = 2` holds
            with both signs for every label `p <= 6`, at every level and for
            all `a, e`.  The words `hat S_p hat S_2 h_1^a` and
            `h_(p-1) hat S_2 h_1^a` are the `b = 1` cases of those theorems.
        - FM-CHK55 (luna_max_eris, fresh code): ACCEPT all six items.
          - Lemma 1: normalization `4^(N+1)/pi^2` and the fold.
          - Lemma 2: the IBP identity (5) coefficient by coefficient;
            Prop. 3.
          - The six polynomial identities and Bernstein minima; the finite
            box is exactly the complement (28,946 values, 10 zeros).
          - Prop. 7 endpoints and the label-7 ray values.
          - The `[U_n] g` identities (1,152 direct cases).  A fresh
            U-positivity screen with `k + a + 2b <= 40` (5,950 profiles)
            found no negative coefficient.
      - **FM-MECH47 (astra_max_ceres; `fm39/mech47_sector123_repro.py`, rerun
        exactly): the {1,2,3} sector is PROVED at every level.  Every FM3
        list with labels `<= 3` satisfies FM3.  In consumer form,
        `phi_r(h_2^al h_1^a hat S_2^b hat S_3^ga) >= 0` for all `r >= 1`,
        `a, b, ga >= 0` and `al <= 2r`.**
        - Theorem 1: if `A, E` are even and `A + E >= al + ga`, then
          `M = E[s^A d^E Z^b (Z+P)^al (Z-P)^ga] >= 0`.  FM3 words satisfy
          this, since `al <= E` and `ga <= A`.
        - Ray form.  With `x = 2 sqrt z`, `y = 2q sqrt z`, one has
          `Z +- P = lam_+- z - 2`, `lam_+- = 4(1 + q^2 +- q) >= 3`.
          Normalize each factor by its value at `z = 1`
          (`f = (beta z - 2)/B`, `h_+- = (lam_+- z - 2)/(lam_+- - 2)`),
          and put `L = al + ga`, `N = (A+E)/2`.
        - Lemma 3 (exponential suppression).
          - On `z >= 2/3` every normalized factor is `>= 3z - 2 >= 0`, and
            AM-GM gives `z >= (3z-2)^(1/3)`.  So the positive part is
            `>= 1/[9(K+b)(K+b+1)]` with `K = L + N/3 + 1`.
          - On `z <= 2/3`, `3z - 2 <= h <= z`, so `|h| <= 2` and
            `sqrt z |h| <= theta = 4 sqrt 2/9 < 1` (maximum at `z = 2/9`).
          - Since `N >= L/2`, each label-3 factor carries a `sqrt z`.  So
            the negative part is at most
            `4 theta^(L-2) (2/3)^(N-L/2) [1/(4(b+1)(b+2)) + (32/225)(5/9)^b]`.
          - A second bound, (5), handles small `N` uniformly in `L`.
        - Lemma 4: every ray integral is `>= 0` once `L >= 22`, `N >= 20` or
          `b >= 42` (exact rational checks and monotonicity).
        - Lemma 5: `L <= 1` is FM-MECH41/44.  `L = 2`, `N <= 1` follows from
          character identities, e.g.
          `Z^2 - P^2 = 1 + U_4(x) + U_4(y) + U_2(x) U_2(y)`.
        - The remaining box is `2 <= L <= 21`, `max(2, ceil(L/2)) <= N <= 19`,
          `b <= 41`.  All 1,848,168 exact values are positive; the least is
          `M(2,2,0,1,1) = 2`.  Evaluation uses an IBP recursion for
          `E[s^(2m) d^(2r) P^j Z^b]` (1,242 bridges) and 528 direct
          consumer checks.
        - Main-agent checks.
          - Every analytic inequality re-derived by hand: `h >= 3z - 2`,
            `f >= 2z - 1`, `w >= 1 - z` on `[2/3,1]`; the `theta`
            maximum; `|f| <= 1 - 2z` on `[0,2/5]` and `|f| <= 5/9` on
            `[2/5,2/3]` (using `beta in [4,8]`); the constant 481/50; the
            monotonicity steps in `N`, `L` and `b`.
          - `fm39/mech47_box_independent_check.py`: 400 random box entries
            agree with an independent `(u,v)`-expansion evaluator.
          - `fm39/mech47_ray_bounds_numeric.py`: bounds (4) and (5) hold on
            3,000 random rays (numerical quadrature).
        - Not claimed: the unconstrained form with `A + E < al + ga` and
          `al + ga >= 3` (not needed by FM3).  `al + ga <= 2` holds without
          the condition.
        - FM-CHK56 (luna_max_eris, fresh code): ACCEPT all five items.
          - The ray normalization.
          - Every pointwise inequality; bounds (4) and (5) at 80-digit
            precision on 48 rays, minimum slacks 0.00028 and 0.00103.
          - The constants and monotonicity steps of Lemma 4.
          - The character identities; 25,000 random box entries with an
            evaluator independent of recursion (10), all positive; (10)
            itself on 1,225 tuples.
          - The consumer translation `A + E - (al + ga) = a + 2r - al >= 0`.
        - FM-SEC128 (luna_max_mercury; `fm39/sec128_ratio_induction_repro.py`,
          rerun exactly): an induction on `al + ga` for the unconstrained
          `b = 0` form.
          - `al + ga <= 2` proved by explicit numerators, e.g.
            `R_(1,1) = 4Q/D_2` with
            `Q = p^4 + 6p^3 + 11p^2 + 6pq^2 + 12p + 4q^4 + 2q^2 + 18`,
            `p = m+k`, `q = m-k`.
          - The unnormalized table is reverse-TP2 on all 35,960 squares.
            That alone cannot drive the recurrence: the constant array has
            zero minors and a negative update.
          - KILL: the ratio-monotonicity rules, both the coordinatewise
            rule and the one oriented by the sign of `al - ga` (witness
            `(m,k,al,ga) = (1,0,0,1)`); normalized TP2 fails (11,120
            negative minors).  Knob: the direction of the ratio bounds.
          - Superseded for FM3 by FM-MECH47; the unconstrained form with
            `al + ga >= 3` stays open and is not needed.
      - **FM-MECH48 (astra_max_ceres; `fm39/mech48_labels4_repro.py`, rerun
        exactly with 32 threads: PASS): the labels `<= 4` sector is PROVED
        at every level.  `phi_r(h_2^al h_3^m h_1^a hat S_2^b hat S_3^ga
        hat S_4^n) >= 0` whenever `al + m <= 2r`, with `b`, `ga`, `n`
        unrestricted.  So FM3 holds for every list with labels `<= 4`.**
        - Uniform contraction lemma (every label).  With `t = sqrt z`,
          `T = [[t, i sqrt(1-t^2)], [i sqrt(1-t^2), t]]`,
          `Q = diag(e^(i phi), e^(-i phi))`, `q = cos phi`, and
          `A_+- = I +- Re R_n(Q) >= 0` on `V_n`, one has
          `tr(R_n(T) A_+-) = U_n(2t) +- U_n(2tq)` and
          `tr A_+- = n + 1 +- U_n(2q)`.  So every normalized Chebyshev
          pair has absolute value `<= 1` on the ray.  The diagonal
          entries of `R(T)` give the strict margin
          `theta_n(rho) = 1 - (1-rho)/(2 C_n)`,
          `C_n = (n+1)(2n^2+4n+3)/3`.
          - Main-agent check of the lemma: `TQ` and `TQ^(-1)` have trace
            `2tq`, and `|tr(UA)| <= tr A` for `U` unitary and `A >= 0`.
        - Schema (every fixed `k`).  Extract `d`, `ds`, `s` or nothing per
          label type (radial resources 1/2, 1, 1/2, 0).  Then
          `delta = 1/(2k^2)`, `rho = 1 - delta`, and Markov's inequality
          gives `W >= 1 - 2k^2(1-z)` near the edge.  Suppression
          (FM-MECH47 Lemma 3 with these constants) leaves an explicit
          finite box `B_k`, and FM3 for labels `<= k` is equivalent to
          positivity on `B_k`.  The constants give
          `1 - theta = Theta(k^-5)` and `H_0 = O(k^5 log k)`.
        - `k = 4` with sharper constants (`rho = 7/9`; `Z +- P`, `Z - 1`,
          `hat S_4` bounds 7/10, 5/8, 2/3): all `H >= 37` nonbase factors
          are covered.  The box has 37,276 records (`N <= 33`, `b <= 77`,
          `H <= 36`): 54,387,664 exact values, one zero, least positive
          value 2.  Also 2,240 shifted-moment bridges and 200 direct
          Catalan checks.
          - Main-agent check of the code: each word is an exact mpz
            combination of table entries, with in-place transforms for
            `Z - 1` and `hat S_4 = Z^2 + Z - 2P^2`.  Every value is
            sign-checked.
          - `fm39/mech48_random_words_check.py`: 400 random labels `<= 4`
            words, evaluated independently, are all `>= 0`.
        - Label 5 (`fm39/mech48_label5_rays_repro.py`, rerun: PASS).  Exact
          ray certificates at `rho = 6/7` prove the labels `<= 5` sector
          for all words with `H >= 71` nonbase factors.  The box below
          that has 1,088,430 core count patterns and is not run.
        - FM-SEC134/135 (luna_max_mercury; `fm39/sec135_labels5_constructor_repro.py`):
          the labels `<= 5` box at `rho = 6/7` (`delta = 1/7`,
          `tau = 17/21`, `J = 352/1225`).
          - Per-type constants certified by 16x16 rational Bernstein
            subdivisions, e.g. `|hat S_4 W| <= 85/147` on `[0, 6/7]`.
          - The many-factor tail is covered for `H >= 68`.
          - The box: 916,895 core profiles, 1,738,408 records,
            227,336,512,420 entries (`N <= 71`, `b <= 1006`).
          - A 32-thread calibration checked 406,203,451 entries (0.18% of
            the box) in 6,838 s: no negative value, and 200 direct
            Catalan matches.  The projected full run is about 44 days.
          - Decision (main agent): the full box is NOT run.  Per-label-bound
            boxes cannot finish the cone, and this one is impractical.
            The labels `<= 5` sector should come from the layer theorems
            (H = 1, H = 2, ...) plus the many-factor theorems.  The
            constructor rerun reproduced the counts exactly (916,895 /
            1,738,408 / 227,336,512,420) and its certificate checks: PASS.
        - Edge obstruction.  With `z = 1 - sigma/M^2`,
          `q = 1 - tau/(2M^2)`, `M = n + 1`, the normalized pair tends to
          `E_+-(sigma,tau) = [F(sigma) +- F(sigma+tau)]/[1 +- F(tau)]`,
          `F(t) = sin(sqrt t)/sqrt t`.  This is an average over the
          uniform sphere in `R^3`, so `|E_+-| <= 1`, strictly at fixed
          `sigma > 0`.  A label-independent quantitative margin, e.g.
          `sup{|E_+-| : sigma >= 1/2} <= 99/100` and its finite-`n`
          version, would make `H_0` logarithmic in `k`.
        - FM-CHK57 (luna_max_eris, fresh code): ACCEPT items 1, 2, 3, 5.
          - The trace lemma and margin (4,999 numeric trace comparisons).
          - The `k = 4` constants and the `N = 1` identities.
          - The consumer translation.
          - The label-5 Bernstein controls.
        - REPAIR (wording) of item 4: the box is a covering superset of
          the analytically uncovered region, not its exact complement.
          2,069 records include 84,415 entries the bound already covers,
          e.g. record `(N,l,m,n,B) = (13,0,0,1,4)` at `b = 3`.  The proof is
          unaffected, since every box entry is evaluated.  The checker
          reproduced the record and entry counts and all cutoffs for
          `H <= 6` plus 600 sampled records.  A fresh 20,000-entry direct
          sample is positive (least value 32).  The full box was rerun by
          the main agent: PASS.
        - FM-SEC130 (luna_max_mercury; `fm39/sec130_labels45_screen_repro.py`,
          rerun exactly in 2 min): independent screens.
          - Labels `<= 4` to total degree 40 (861,016 words) and labels
            `<= 5` to degree 30 (295,668 words): no negative value.
          - The factor catalogue through label 12 was verified, e.g.
            `+5 = s(Z^2 - ZP - P^2 + 2P - 1)`; 512 Catalan bridges.
          - Tightness versus the largest label `n`.
            - One label `n`, the rest in {1,2,3}: the normalized minimum
              is `2/mu(0, ceil(n/2))`, from the degree-boundary word
              `(-1)^n (-+n)`, whose value is exactly 2 (distance 0, inside
              H_AC_q).
            - Two labels `n`: the minimum stays in `[2/5, 2/3]` through
              `n = 12`.
      - **FM-MECH49 (astra_max_ceres; `fm39/mech49_uniform_margin_repro.py`,
        rerun exactly): a label-uniform many-factor theorem.**
        - Theorem 1 (edge margin, every `n >= 3`, `M = n + 1`):
          `|U_n(2 sqrt z) +- U_n(2q sqrt z)| / [M +- U_n(2q)]
          <= 1 - min{M^2(1-z), 1/2}/125`.  In particular
          `c = 249/250` on `z <= 1 - 1/(2M^2)`.
          - Proof: the FM-MECH48 trace identity; middle weights of
            `Sym^n` keep a fixed fraction of the mass of `A_+-` (rotation
            energy of order `n^2`); scalar sine bounds for `X = M sin(alpha)
            >= 4`.  No limiting approximation is used.
          - Main-agent check: 120,000 random high-precision tests,
            concentrated near the edge, with no violation.
        - Lemma 4 (positive region): Markov's inequality gives
          `W_n >= 1 - C M^2 (1-z)` near `z = 1`.
        - Theorem 5.  Every consumer word with `H` factors of label `>= 3`
          (signs free), at every level and with any number of 1's and 2's
          (all `N`, `b`), satisfies FM3 once a weighted count `T` of those
          factors reaches `T_0(k) = O(log k)`, `k` the largest label.  In
          plain count `H >= H_0(k) = O(k^2 log k)`, e.g. `H_0(5) = 19,612`
          and `H_0(10) = 78,848`.  If every label is `>= eps (k+1)`, then
          `H >= T_0/eps^2` suffices.  (At `k = 5` the FM-MECH48 cutoff
          `H >= 71` is stronger.)
        - (Residual-filter wording repaired by ADV-3: apply LS, LR4 and
          LLm-S to the whole non-fundamental word, and keep the original
          list.)
        - FM-CHK59 (luna_max_eris, fresh code): ACCEPT items 1-3.
          - Theorem 1: all scalar estimates, the sine certificate (slack
            3523329/167772160) and the rotation-energy lemma.  High-precision
            checks: 6,000 random signed comparisons (1,400 edge-biased) and
            1,788 removable endpoint limits, minimum sampled slack
            `9.2e-122` (at the edge, where both sides tend to 1).
          - Lemma 4.
          - Theorem 5: the weighted count, the monotonicity in `N` and
            `b`, and the recomputed `T_0`/`H_0` table.
        - REPAIR item 4: Prop. 7's factor-count front predates Theorem Q2+
          (FM-MECH50) and Theorem B2 (FM-MECH51).  `H = 1, b = 1` is now
          covered, e.g. `phi_2(hat S_7 hat S_2 h_1^7) = 88`; the next open
          slice is `b >= 3` (FM-MECH52).  The `(3,5)` maximum-label front
          stands.
        - Proposition 7 (residual of the full cone).  The residual is the
          union over `k` of explicit finite sets `B_(k,new)`, minus the
          proved strata.  Necessary conditions for a residual word: list
          distance `d >= 3` with some label below `d`; LS, LR4 and LLm all
          fail; no OL or proved (E) region applies.
          - First front by factor count: one label `n >= 7` with one
            `hat S_2` and a fundamental suffix,
            `phi_r(hat S_n hat S_2 h_1^a)` and `phi_r(h_(n-1) hat S_2 h_1^a)`
            at `d >= 3`.  This is exactly the q = 2 plus-W sign of (E)
            (FM-SEC127..133).
          - First front by maximum label: `k = 5`, `2 <= H <= 70`, first
            core pattern `(3,5)`.
      - **FM-MECH50 (astra_max_ceres; `fm39/mech50_q2plus_repro.py`, rerun
        exactly in 36 s): Theorem Q2+, the q = 2 plus-W sign of (E), PROVED
        uniformly.  For every `r >= 1`, `a >= 0`, `p >= 3`:
        `phi_r(hat S_p hat S_2 h_1^a) >= 0` and
        `phi_r(h_(p-1) hat S_2 h_1^a) >= 0`.  This closes the first front of
        the FM-MECH49 residual and the FM-SEC35/38/44/127..133 item.**
        - Residual statement: `min(a,e) >= 3`, list distance `n = N - j >= 3`,
          `kappa = 2j - N >= 5` (label `p = kappa + 2 >= 7`).  The other cases
          are the strips, H_AC_q (`d <= 2`) and FM-MECH45 (`p <= 6`).
        - The slack times `H` is the quadratic form
          `A c_j^2 - Delta beta c_j c_(j+1) + gamma c_(j+1)^2`, with
          `A > 0` explicit (as in FM-SEC129).
        - Coverage lemma.  The discriminant `J(xi)`, `xi = Delta^2`, is a
          cubic with `[xi^3] = (kappa+3)^2 > 0`, positive `[xi^2]` and
          `J(0) < 0`, hence convex on `xi >= 0`.
          - For `n <= 2(kappa+1)^2`: `-J(4O)/4` has an explicit
            positive-coefficient form, so `J < 0` below the outer
            threshold `Delta^2 = 4O`.
          - For `n >= 2(kappa+1)^2`: `-L^6 J(xi_c)`, with `L = 3(kappa+1)^2`
            and `xi_c = N^2 (L-2)^2/L^2`, has 133 positive coefficients
            after `n = n_0 + u`, `kappa = 5 + v`.  Outside the central
            region, `Delta^2 < xi_c`.
          - So every residual row is PD, outer or central.
        - Outer lemma: the Bernstein coefficients of the form on
          `[0, R_OL]` are certified positive after `z = O R^2/n^2`
          (cubic certificates with 21 to 27 positive terms).
        - Central lemma: the Bernstein coefficients on the OL central ratio
          intervals (even `t`: `[0,1]`; odd `t`: `[l, (kappa+2)/kappa]`) are
          certified positive (140 to 406 terms each).  The odd upper bound
          `(kappa+2)/kappa` is new: it comes from the paired-root product
          of `K_t(x)/x`.
        - Also 819 Catalan bridges against direct consumer values.
        - Main-agent check of the logic: the convexity, chord and region
          arguments are as stated.
        - FM-CHK58 (luna_max_venus, fresh code): ACCEPT all six items.
          - The reduction: 252 direct cases, plus 75 parity zeros.
          - The fold invariance.
          - `A > 0`; the discriminant data, with the full factorization of
            `J(0) < 0` and the five positive `p_h`; the 133-term
            certificate.
          - The OL continued-fraction bound at this index; the
            Bernstein certificates.
          - The new odd-`t` bound `rho < (kappa+2)/kappa`, from Krawtchouk
            duality and the paired-root product (300 central rows checked).
          - An independent screen of all 398,432 residual rows with
            `a, e <= 120`: all positive, minimum slack 64.
        - Next arbitrary-label background (Ceres):
          `E[(x-y)^e (x+y)^a (x^2+y^2-2)^2 U_p(x)] >= 0`, i.e. two
          `hat S_2` factors.
      - **FM-MECH51 (astra_max_ceres; `fm39/mech51_b2_repro.py`, rerun exactly
        in 30 s): Theorem B2 PROVED.  `phi_r(hat S_p h_1^a hat S_2^2) >= 0`
        and `phi_r(h_(p-1) h_1^a hat S_2^2) >= 0` for every `p`, `r`, `a`.
        Also an explicit cutoff, uniform in the label and in `b`.**
        - The `b = 2` slack is a quadratic form in `(c_j, c_(j+1))`
          (degrees 6, 7, 8 in `Delta`), via the half-angle rows
          `4 delta(c) + 4 delta(c^(1)) + delta(c^(2))`.
        - The Q2+ three-region method applies with these changes:
          - the discriminant is a quintic in `Delta^2`; coverage uses its
            six Bernstein coefficients on `[0, 4O]` and on `[0, xi_c]`;
          - the outer region is reparametrized by `x = q - sqrt S`,
            `y = q + sqrt S`, with degree-3 Bernstein coefficients
            (degree 2 has negative entries);
          - the central intervals are those of Q2+, including the odd
            bound `(kappa+2)/kappa`;
          - `t = 0, 1` use explicit Krawtchouk vectors.
          All certificates are fixed polynomials with positive
          coefficients; 299 Catalan bridges.
        - Lemma 4 (every `b`): with `w = z^nu sqrt((1-z)(1-q^2 z))` and
          `f = beta z - 2`, one has
          `int w |f|^b (1-z) <= 4 J_b/(nu + b - 2)`, by a reflection
          pairing and IBP.  So
          `int w f^b W >= [1 - 4 ||W'||/(nu+b-2)] J_b`.
        - Theorem 5 (every `b`).  The FM-MECH49 trace bound gives
          `|z W| <= 1`; Markov's inequality twice gives
          `||W'|| <= 4 (p+1)^4`.  So both families are `>= 0` whenever
          `nu + b >= 16 (p+1)^4 + 2`, with `nu = r + (a + eps)/2`; also
          `nu = 1` by characters.  Main-agent check of the Markov chain:
          correct.
        - Remaining in the `H = 1` layer: `p >= 7`, `nu >= 2`,
          `3 <= b <= 16(p+1)^4 + 1 - nu`, list distance `>= 3`.  This is a
          finite `b`-range for each label, and the missing strength is
          uniformity in `b` there.
        - FM-CHK60 (luna_max_venus, fresh code): ACCEPT items 1-4 and 6.
          - The `b = 2` reduction; the recurrence quadratic and `A > 0`.
          - The quintic-discriminant coverage, all six plus six Bernstein
            checks.
          - The central and outer certificates (17,980 outer rows against
            the continued-fraction interval) and `t = 0, 1`.
          - An exact screen of all 231,624 folded residual rows with
            `a <= 100`: positive, minimum 44.
        - REPAIR item 5 (cutoff): Lemma 4 is stated for `b >= 1`, so
          `b = 0` needs its own line.  Either IBP gives
          `int w (1-z) <= 2 J_0/nu` and the same conclusion, or `b = 0` is
          Theorem OL.  The reflection, IBP and Markov steps are accepted
          (378 quadrature cases, least normalized slack 2.73; 66 exact
          reflection checks).
      - FM-MECH52 (astra_max_ceres; `fm39/mech52_h1_layer_repro.py`, rerun
        exactly): the `H = 1` layer, reduced to one polynomial family.
        - Theorem 1 (PROVED): every `g_(e,a,b)` with `min(a,e) <= 1` is
          U-positive, for every `b` and every label.  Proof: an
          antisymmetric two-variable cone `Psi_(i,j)`, `i > j`, preserved by
          multiplication by `Z = U_2(x) + U_2(y)` and by `x + y` (indices
          reach the diagonal but never cross it).
        - Theorem 3 (PROVED): the cutoff improves to quadratic.  Both
          families are `>= 0` once `nu + b >= K`, with
          - `K = 2p^2 + 2` for `hat S_p`, `p` even;
          - `K = 2p^2 - 2` for `h_(p-1)`, `p` even;
          - `K = 4p^2 - 2` for either family, `p` odd.
          Ingredients: a weighted reflection lemma and the FM-MECH49 trace
          bound with Markov's inequality.  E.g. 194 instead of 65,538 at
          `p = 7`.
        - KILLs (knob: extra growth under insertion, or raywise
          positivity).
          - A `U_2`-transfer bound fails by -6; doubling fails by 424,528.
          - A negative ray at distance 4 (`-17093/495495`, consumer value
            11128), and an unbounded negative-ray family at distance
            `~ p^2` (limit via `sinc(sqrt y)` and density `y e^(-y/20)`).
          So the remaining region needs angular cancellation.
        - Exact distance formula (5): with `T = a + e`, `D = T + 2b`,
          `p = D - 2d`,
          `F_d = sum_(j+l <= d) c_(2j) C(b,l) mu_(j,l)
          [C(D-2j-2l, d-j-l) - C(D-2j-2l, d-j-l-1)]`, where
          `mu_(j,l) = E[y^(2j) (y^2 - 2)^l]`.  It is a polynomial of total
          degree `<= d` in `T`, `(a-e)^2`, `b`, independent of the label.
          Its leading coefficient in `((a-e)^2)^d` is `1/(d!(d+1)!)`.
        - Remaining `H = 1` box: `a, e >= 2`, `b >= 3`, `p >= 7`, and
          `4 <= d <= 4p^2 - (p+1)/2 - 3` (`p` odd) or about `2p^2` (`p`
          even).  The target is `F_d >= 0` there.
      - FM-MECH53 (astra_max_ceres; `fm39/mech53_distance45_repro.py`, rerun
        exactly in 13 s): the `H = 1` layer at distances 4 and 5, every
        label and every `b`.
        - Theorem (PROVED): `F_4, F_5 >= 0` on `a, e >= 2`, `b >= 3`.
          - Large-`T` region (`T >= 7b - 3`, resp. `6b - 1`): a square
            completion plus remainders with positive coefficients.
          - Complementary region: Bernstein certificates in `v` for each
            `[u^k] G_d`.
          - Finite corners: 90 and 557 cases, all positive, least values
            112 and 256.
          - So `H = 1` is closed through distance 5; distance 6 holds for
            `b >= 6`, `T >= (11b-1)/2` (a mixed-term square repair).
        - Proposition 2: `G_d` is monic in `u` with an explicit penultimate
          coefficient for every `d`.  This explains the boundaries near
          `7b`, `6b`, `11b/2`.
        - Proposition 3 (every `d`): the leading homogeneous part `L_d` is
          `>= 0`, via
          `sum_d L_d z^d = exp(Tz) E_y[exp((b - T/2) y^2 z) cosh(y sqrt(uz))]`
          and, at `b = 0`,
          `L_d(T,u,0) = (1/(d+1)!) sum_j T^(d-j) H_j(sqrt u; T)^2/j!` (monic
          Hermite squares).  Also
          `L_d(T,u,b) = sum_h (2b)^h/h! L_(d-h)(T - 2b, u, 0)`.
        - KILLs.
          - The square-remainder template fails at `d = 6`; the
            `T = 11b/2` defect is `-(5005b^2 + 26100b - 1628)/4`.
          - No fixed relative margin `F_d >= delta L_d`: a family below
            the cutoff has `|F_d|/L_d -> 0`.
          Knobs: coefficientwise positivity of a particular remainder,
          and a fixed Gaussian margin.
        - Open: `d >= 6` uniformly.  A per-`d` ladder would not finish.
      - FM-CHK62 (luna_max_eris, fresh code): ACCEPT ADV-1 Lemma 1 (the
        distance-3 identity and certificate), FM-MECH52 (antisymmetric cone,
        weighted reflection bound, quadratic cutoff) and FM-MECH53 (distances
        4, 5; the 90 and 557 corners, least values 112 and 256; the
        distance-6 subregion; the penultimate coefficient; the Hermite-square
        leading part).
      - FM-MECH54 (astra_max_ceres; `fm39/mech54_genfun_obstructions_repro.py`,
        rerun exactly): exact distance generating function; four KILLs.
        - Prop. 1: `F_d = [z^d] (1-z) E_y[(1+z-sqrt z y)^e (1+z+sqrt z y)^a
          (1+z^2+z y^2)^b]` for `d <= D/2` (1,410 Catalan checks).  This is
          the half-angle determinant formula in generating-function form:
          on `|w| = sqrt z`, `1 + z +- sqrt z y = |1 +- w|^2` and
          `1 + z^2 + z y^2 = |1+w^2|^2 + 2|w|^2`.
        - KILLs, all on actual integer consumer profiles:
          - (Prop. 2) The finite background translation has alternating
            Christoffel-weighted terms; the quotient by the shifted `b = 0`
            series is `1 + 6z + 18z^2 - 46z^3` at `(e,a,b) = (3,11,3)`.
            Even `d`-dependent weights over the same shifted basis need a
            `-216` coefficient.
          - (Prop. 3) No nonnegative reduction to `b = 0` profiles that
            holds for all output labels `>= 7`: separating functional
            value `-5608` at `g_(6,8,3)`.
          - (Prop. 4) No diagonal Krawtchouk-square representation with
            imbalance-independent weights (dual value `-2/5` at
            `T = 13`, `b = 3`, `d = 6`).
          - (Prop. 5) Newton positivity in `b` fails (`-32` at
            `T = 14`, `u = 0`, `d = 6`).
          Knobs: simultaneous-output reductions, a single shifted
          background, diagonal squares, Newton growth.  Output-reindexed
          reductions and mixed squares remain untested.
        - Main-agent remark: the min(a,e) <= 1 cone (FM-MECH52 Thm. 1)
          does not extend directly to `e >= 2`.  Multiplying by `d` maps
          the antisymmetric cone to signed symmetric combinations, and
          squares are not U-positive (`(U_1 - U_3)^2` has
          `[U_4] = -1`).
        - Stall: two rounds (FM-MECH53, 54) without a uniform-in-`d` step.
          The demand sheet on this chain is exact (all `d`, `b >= 3`,
          integer `a, e >= 2`, both families), so nothing restates.
          Effort is split: an advisor on the `H = 1` strategy, and Ceres
          on the `H = 2` layer ((E) for `q >= 3`).
      - FM-MECH55 (astra_max_ceres; `fm39/mech55_E_energy_repro.py`, rerun
        exactly): (E) in two new regions, uniform in the gap `q`.
        - Theorem: with `sigma = N + 2`, `Delta = a - e`, `X_k = 2k - N`,
          `s = i - j >= 4`, `m = floor(s/2)` and `4|Delta| <= sigma`, both
          signs of (E) hold, at every level and every `q`, in two regions:
          - the band `sigma/4 <= X_j <= sigma/2`;
          - the long interval `0 <= X_j <= sigma/2` with
            `3m(X_j + 2m - 2) >= 16 sigma`.
          Both signed slacks are positive definite there.
        - Tools.
          - A fixed quadratic energy `H` on consecutive coefficients, with
            an exact square decrement on the right half.
          - A rotation-then-contraction factorization of the recurrence
            step, with two-step operator-norm bounds.
          - An endpoint comparison.
          - One fixed Bernstein certificate (bidegree (9,5), 60
            coefficients, minimum 462163/32).
          - 39,246 energy checks on 1,891 rows; the regions contain 34,034
            row pairs, least slack 36; 705 Catalan reduction pairs.
        - KILL (knob: separate optimization of the two endpoint directions):
          the endpoint envelope fails at `(a,e,j,i) = (48,8,29,33)`
          (squared deficit < 0), while the actual slack is positive.
        - Open: coverage of the complement (region (10) of the report)
          intersected with the failure sets of the short-arc, RF, LD and
          metric criteria.
      - FM-MECH56 (astra_max_ceres; `fm39/mech56_E_propagation_repro.py`,
        rerun exactly): (E) reduced to an explicit minimal-counterexample
        set.
        - Prop. 1 (direction-dependent propagation): if the edge
          `psi_k -> psi_(k+1)` crosses the ray opposite `psi_j`, an exact
          identity writes the parent slack as a weighted sum of the two
          child slacks plus a correction `G_jki`.  So (E) on both shorter
          intervals together with `G_jki >= 0` gives (E) on `[j,i]`.  This
          is an induction on the gap, with OL, EQ1 and Q2+ as base cases.
        - Prop. 2: an opposite vertex `psi_k = -lambda psi_j` cannot occur in
          a minimal-gap counterexample.  The centre `j = N/2` is direct.
        - Prop. 3: with FM-MECH55's energy and the actual drop `D_j - D_i`,
          a scalar `J_ji >= 0` implies (E).
        - KILLs (knob: replacing child margins by zero):
          - `G_jki < 0` occurs on actual rows, e.g. `(31,19,26,28,39)`;
          - the drop condition without `D_i` fails at a four-step crossing
            `(172,45,110,114)`.
        - Theorem 5: (E) is equivalent to (E) on an explicit set `R`
          (folded, outside FM-MECH55, arc `>= pi`, no opposite vertex, every
          crossing with `G < 0`, and `J < 0`).
        - Census (main-agent rerun): 729 rows, 117,471 long-arc pairs, all
          covered.  Gluing covers 115,682, opposite vertices 2,005, and the
          `J_ji >= 0` criterion alone all 117,471.  `R` is empty in this
          range.  Next target: prove `J_ji >= 0` on long arcs uniformly.
      - FM-MECH57 (astra_max_ceres; `fm39/mech57_J_propagation_repro.py`;
        default part and the `--large --threads 32` census both rerun exactly by
        the main agent, about 2 minutes): the joint-energy criterion
        `J_ji >= 0` on folded long arcs.
        - Exact screens (folded `a >= e+2`, `e >= 3`, `i - j >= 4`) with
          `R = ` normalized `J`, where `J >= 0` iff `R >= 1`.  No failure in
          any of them:
          - `a, e <= 200`: 104,658,739 long pairs, minimum `R = 2.0546`;
          - 256 random rows up to 1000: 52,916,472 pairs, minimum `2.0694`;
          - lopsided rows up to `a = 4000`: 1,654,937 pairs, minimum `1.9454`.
          All pairs also satisfy (E).
        - Proposition 2 (PROVED, uniform outer propagation): past the turning
          point `X_i^2 >= C`, `J(j,i) >= 0` implies `J(j,h) >= 0` for every
          later endpoint `h`.  Proof: a continued-fraction bound from the
          boundary `c_(n+1) = 0`, the exact square decrement, and OL.
        - Proposition 3 (KILL of the unfolded version; knob: folding): for
          every `m >= 28` an explicit gap-4 long arc has `J < 0` before
          folding, with `R -> 0`; after folding it is a short arc, already
          covered.  Certificate: degree 16, all 17 coefficients positive.
        - Corollary: for each folded row and anchor, the remaining work is
          `J` at the endpoints from the first long crossing `i_0` up to the
          turning point.  The forward increment was nonnegative in all
          54,347,507 tested steps, but this is not proved.
      - FM-CHK63 (luna_max_pluto, fresh code): ACCEPT FM-MECH55 (energy,
        contraction bounds, band and long-interval regions), FM-MECH56
        (propagation identity, opposite-vertex elimination, `J` criterion,
        reduction to `R`) and FM-MECH57 Props. 2, 3 (outer propagation, the
        unfolded counterfamily).  The kill witnesses were reproduced.
        Wording note: "necessary" in the FM-MECH57 corollary is relative to
        that propagation scheme.
      - FM-MECH58 (astra_max_ceres; `fm39/mech58_J_central_repro.py`, rerun
        exactly in 51 s): central-anchor coverage of `J`.
        - Theorem 4 (PROVED): for every folded row with
          `4|a - e| <= N + 2` and `2X_j <= N + 2`, `J(j,i) >= 0` at EVERY
          endpoint `i >= j + 4`.  So both signs of (E) hold there, uniformly
          in the gap `q`.  In two-label coordinates `X_j = p - q`.
        - Ingredients:
          - Lemma 2: two-step contraction bounds, uniform in the gap, for
            `i - j >= 6` (with `X_i <= sigma/2`) and `i - j >= 8` (with
            `X_j <= sigma/2`); fixed Bernstein certificates of degree
            (6,3), 28 coefficients each.
          - Lemma 3: gaps 4 and 5 by exact binary quadratic forms (diagonal
            and determinant certificates with 660 to 3,750 positive
            coefficients).
          - The FM-MECH55 band.
          - A fixed finite remainder of 475 cases (`N <= 53`), least ratio
            182329/26550.
          - Rerun: 100,867 region pairs, 48,573 of them beyond FM-MECH55.
        - Open (Prop. 5): `J` on the folded long arcs with `4|a - e| > N + 2`
          (unbalanced rows) or `2X_j > N + 2`, from the first long endpoint
          up to the turning point.  The endpoint-envelope obstruction
          `(48,8,29,33)` (actual `J >= 0`) shows the endpoint directions
          must stay correlated.
      - FM-MECH59 (astra_max_ceres; `fm39/mech59_E_far_anchor_repro.py`,
        rerun exactly in 46 s): the far-anchor branch of (E), CLOSED.
        - Lemma 1 (a square decrement for `D`, Sonin-type): with an explicit
          rate `rho(X_k)`, `D_(k+1) <= rho(X_k) D_k` on the right half.  This
          is an exact identity whose remainder is a square, since
          `(sigma - omega)(sigma + omega) = Delta^2`.
        - Lemma 2: the decay rate combined with Theorem LD's chord bounds
          gives (E) for every far-anchor interval of length `>= 5`, and for
          length 4 when `5X_j >= 3 sigma`.
        - Lemma 3: the remaining gap-4 region has PSD binary forms, by fixed
          certificates (330 to 7,605 positive coefficients).
        - Theorem 4 (PROVED): (E), both signs, for every far anchor
          `2X_j > N + 2`, at every endpoint and every gap.  The other cases
          use Theorem RF's outer region and 100 exact cases (least slack
          719,376); the endpoint `i = N + 1` uses Newton's inequality.
          Census: 170,240 pairs, all covered.
        - Remaining for the two-label layer: unbalanced rows at central
          anchors, `8e < 3N - 2`.
        - KILL (knob: bounding the two chords separately): at
          `(4000,4,2016,2051)` the LD bound exceeds the drop by a factor
          above `sqrt 10`, while `J >= 0`.
      - FM-SEC138 (luna_max_mercury; `fm39/sec138_coverage_atlas_repro.py`;
        phase 0 rerun exactly in about a minute, and the `--match` interface
        spot-checked): an exact coverage atlas of the proved strata.
        - The matcher keeps the original signed list and returns every
          matching theorem tag.  It is conservative: some proved criteria are
          not encoded (parts of the (E) regions such as root-free, binary
          form and LD/metric chords; the full EQ1 row; G0E4; the finite
          H-only box).
        - Phase 0 (labels and sum `<= 30`): 1,114,614 lists, 129,186
          unmatched.  Phase 1 (at most 3 labels `>= 3`, labels `<= 40`, sum
          `<= 60`): 78,037,215 lists, 10,189,015 unmatched.  Every unmatched
          list has a positive exact value; the minimum is 32, at
          `(+1,-2,+2,-3,-5,-5)`.
        - The smallest unmatched lists have total label sum 18, always a
          label 5 together with 3 or 4 (two to three core factors), at
          distance 3 or 4, e.g. `(-1^8, +2, -3, -5)` (value 11014).  So the
          labels `<= 5` front is the first real residual.
        - Theorem LLm has the most matches, then LLm-S, H_AC_q (d <= 2) and
          labels `<= 4`.
        - Main-agent note: the FM-SEC135 labels `<= 5` box (227.3B entries)
          was projected at 44 days with that evaluator.  The FM-MECH48
          incremental C++ evaluator runs about 1,000 times faster (54.4M
          entries in 0.6 s on 32 threads), which would put the box at
          roughly 40 minutes.  Delegated as FM-SEC139 (a port, run by the
          main agent with checkpoints).
      - FM-MECH60 (astra_max_ceres; `fm39/mech60_E_unbalanced_repro.py`,
        rerun exactly in 45 s): unbalanced rows at central anchors.
        - Lemma 1: a joint local metric `Q_X` with `D_k = Q_(X_k)(psi_k)` and
          `Q_X - Q_Y >= 0`, giving a smaller-root criterion for both signs
          of (E).  It keeps `W` as one determinant.
        - Lemma 2: the binomial-decay energy
          `E_k = u_k^2 + V_(X_k) v_k^2 = D_k + 2(sigma - X_k) v_k^2`,
          with `P_X = C - X^2` and `V_X = C + 2 sigma - X(X+2)`, decays by
          exact factors `beta_(k+1)/beta_k`.  This gives the parameter-only
          criterion
          (BE) `(omega + Y) V_X beta_i <= (omega - Y) P_X beta_j`, `Y < omega`.
        - Theorem 3 (PROVED): both signs of (E), at every endpoint including
          past the turning point, in two regions:
          - C1: `X_j <= omega/4` and `X_i >= omega/2`;
          - C2: `X_j <= 2 omega/5` and `X_i >= 5 omega/8`.
          The proof uses checkpoint constants, the Sonin identity and the
          fixed energy.  The fixed remainder `omega < 64` has 490,037 pairs
          on 1,494 rows, least slack 67.
        - Audit: (BE) covers all 998,772 inner long arcs; the earlier Sonin
          parameter bound misses 7 of them, all covered by (BE).
        - Remaining for the two-label layer, three ranges:
          - `0 < X <= omega/4` with `Y < omega/2`;
          - `omega/4 < X <= 2 omega/5` with `Y < 5 omega/8`;
          - `X > 2 omega/5`.
          Plus `J` at the first turning endpoint outside C1, C2.
      - FM-CHK64 (luna_max_venus, fresh code).
        - ACCEPT FM-MECH58: all certificates regenerated (K6, K8, the band,
          gap 4/5).
        - ACCEPT FM-MECH60 for the regions C1, C2.
        - REPAIR FM-MECH59: Lemmas 1-3, the case division, the RF step and
          the 100-case remainder are accepted.  But the terminal case
          `i = N + 1` cited Newton's inequality, which does not apply
          directly to a mixed-sign row.  Replacement (verified by the main
          agent symbolically and on all integer rows with `a < 40`): with
          `b = a - e >= 2`,
          `4(N-2) D_3 - (N+1) c_3^2
          = e(b+e)(b^4 - b^2 + 12b(e-1) + 12(e-1)^2)/3 > 0`.
          By `c_(N-k) = (-1)^e c_k`, `D_(N+1) = 0` and integrality, this
          gives `D_j >= |c_j|` at `j = N - 3`.  So the far-anchor theorem
          stands with this step replaced.
      - FM-MECH61 (astra_max_ceres; `fm39/mech61_E_checkpoint_repro.py`,
        rerun exactly in 10 s).
        - KILL of (BE) (knob: a parameter-only strengthening): for every
          `h >= 4` an explicit long arc in R3 makes (BE) fail, and the
          smaller-root condition also fails (at `h = 8`).  `J > 0` there.
        - Theorem 2 (PROVED): `J(j,i) >= 0` at every endpoint, including
          past the turning point, in two checkpoint regions, one for `e = 3`
          (T3) and one for `e >= 4` with `X <= omega/2` (T4).  This covers
          the whole infinite (BE)-failure family.  For `e = 3` with
          `X^2 >= 3 sigma - 2`, RF gives (E) directly, so the turning
          obligation is closed at `e = 3`.  Fixed remainder `omega < 64`:
          1,105,010 pairs, least slack 67; 133,792 (T3) and 317,222 (T4)
          instances.
        - Lemma 3: on the inner endpoint lattice, the (BE)-successful
          endpoints form one consecutive interval (a unimodality argument
          via a cubic).  So R1 and R2 reduce to (BE) at the first long
          endpoint `i_0`.
        - Remaining for the two-label layer (b = 0):
          - (BE), or a row-based bound, at the initial crossing `i_0`
            before the checkpoint;
          - the larger R3 anchors with `e >= 4`, `omega >= 64`,
            `omega/2 < X <= sigma/2` and `X^2 < 4(e-1)(a+4)`.
      - **FM-MECH64 (astra_max_minerva; `fm39/mech64_two_core_distance5_repro.py`,
        rerun exactly in 6 s): two core labels on every {1,2} background,
        through distance 5.**
        - Combined half-angle identity: with `P_h = (1+z^2)^h P` and
          `j = (N+n-m)/2`, `i = j + m + 1`,
          `R_(n,m) = sum_h C(b,h) 2^(b-h) (D_(j+h)(P_h) - D_(i+h)(P_h))` and
          `c_(n,m) = sum_h C(b,h) 2^(b-h) W_(i+h,j+h)(P_h)`.  At fixed list
          distance `delta` these become polynomial in `b`, with
          `R = T_delta - T_(delta-m-1)`.
        - Theorem (PROVED): both compatible signs for two labels `n, m >= 3`,
          every `b >= 1`, every level and all labels, at distances 3, 4, 5.
          - Diagonal expansions `T_delta = sum_k A_k^(delta) q_k^2` in formal
            Krawtchouk polynomials with `M = N - 2b`, e.g.
            `A_delta = 1/(delta+1)` and
            `A_(delta-1) = (N + 2 delta b - 2 delta + 2)/(delta(delta+1))`.
          - Mixed square completions such as (9) at distance 5.
          - 18 coefficient certificates; `b`-tails (12, 224, 1,944
            polynomials); finite corners (376 and 13,368, minimum slacks 14
            and 31); 16,746 character bridges.
        - Theorem (PROVED): the strip `min(a,e) = 0`, every distance, every
          `b`, all labels, via the antisymmetric cone (`s` and `Z` preserve
          it).
        - With H_AC_q (`d <= 2`), the two-core layer is closed through
          distance 5 on every {1,2} background.  This removes the smallest
          uncovered lists of the FM-SEC138 atlas (distances 3, 4).
        - KILLs.
          - Rank-one summand positivity: `(R,W) = (-1,0)` for the `P_1` part
            of `P = (1-z^2)^2`, `n = m = 3`, against the total `(7,0)`.
          - Insertion on the whole two-label moment cone: the abstract
            `G = 1 + U_3 U_3 - U_2 U_2` has every two-label inequality, but
            `ZG` has `R_(3,3) = 0 < 2 = c_(3,3)`.
      - **FM-MECH63 (astra_max_vulcan; `fm39/mech63_one_label_small_background_repro.py`,
        rerun exactly): one arbitrary label on ANY background of labels
        `<= 4`, above a quadratic cutoff.**  For a word with one label `n >= 5`
        (either sign) and any multiplicities and signs of the labels 1..4,
        put `V = (c_1^+ + c_1^-)/2 + c_2^-`, `b = c_2^+`,
        `J = c_3^+ + c_3^- + c_4^+ + c_4^-` and `B = V + b + J`.  Then FM3
        holds whenever `B >= max(384(n+1)^2, 2^17)`.  Every normalized ray
        integral is positive.
        - Pair estimates on the ray, with `u = 1 - z`:
          - `|R_3|, |R_4| <= 1 - u/3`;
          - `R_3, R_4 >= 1 - 8u` for `u <= 1/4`;
          - `|R_n| <= 4 sqrt z` for odd `n`;
          - `R_n >= 1 - 2(n+1)^2 u`, by Markov's inequality.
          Main-agent check: 160,000 high-precision samples, no violation.
        - Key step: the consumer's parity constraints supply one factor `z`
          uniformly in `n` for the whole product:
          `z^V |R_n prod R_small| <= 12 z exp(-Cu/3)`, `C = V + J`.
        - What it gives: one huge label against arbitrarily many 3s and 4s
          is covered without the weighted many-factor threshold.  The
          remainder is a finite set for each `n` (`B < 384(n+1)^2`).
          Independent check pending (FM-CHK65).
      - FM-MECH66 (astra_max_ceres; `fm39/mech66_E_e3_anchor_repro.py`, rerun
        exactly in 8 s).
        - Theorem 1 (PROVED): (E), all signs, for every row with
          `min(a,e) <= 3`.  The `e = 3` branch is closed, using the exact RF
          formula, (BE) before the T3 checkpoint, and T3.
        - Lemma 2 (arc speed): a diagonal contraction turns a vector by a
          bounded angle, so a long arc forces `kappa(y - x) > 7/5`, i.e.
          enough binomial decay.  Constants come from rational Taylor
          bounds.
        - Theorem 3 (PROVED): both signs at every endpoint when `kappa >= 15`
          and the anchor `X <= 3 omega/5` (outside the very central part),
          including through the turning point.  An exact finite part of
          seven `e = 7` rows (15,341 pairs, minimum slack about `1.38e13`).
        - Remaining for (E) with `q >= 3`:
          - `e in {4,5,6}` inside the FM-MECH61 residual;
          - very central anchors (`kappa >= 15`, `omega X < 6(sigma+1)`,
            `Y < omega/2`);
          - larger anchors (`kappa >= 15`, `3 omega/5 < X <= sigma/2`,
            `X^2 < 4(e-1)(a+4)`).
      - **FM-MECH68 (astra_max_juno; `fm39/mech68_distance34_all_lists_repro.py`,
        rerun exactly in 33 s): FM3 holds for EVERY list at list distance
        `delta <= 4`, whatever its labels, signs, background and number of
        factors.**  (`delta <= 2` is FM-MECH38.  This settles FM-MECH39's open
        item "mixed d = 3 with two or more labels 1".)
        - Prop. 1 (count reduction): with the maximum label distinguished,
          `EVEN = 2 F_delta`,
          `F_delta = sum_S eps_S m(mu_S) K_(delta - w(S)/2)(mu_(S^c))`,
          with `K_j(nu) = [v^j](1-v) prod_i (1 + ... + v^(nu_i))`.
          Non-distinguished labels `> delta` clip to `delta + 1`, and each
          `-2` can be replaced by `+1, -1`.  So `F_delta` is a polynomial in
          the signed counts `t, x` (labels 1, 2), `b` (`+2`), `k` (labels
          `>= 3`), `c, y` (label 3), `f, z` (label 4).
        - Prop. 2 (`delta = 3`): `F_3` equals
          `P_3^2/4 + C_2 P_2^2 + C_1 x^2 + C_0 + y P_3 + (y^2 - c)/2`, with
          explicit `C_i` and shifted Krawtchouk `P_j` (`M = t - 2b`).  A square
          completion gives a lower bound whose coefficients are
          positive-coefficient polynomials after `t -> t+11`, `b -> b+6` or
          `k -> k+4`.  The remaining box has 7,869 count profiles, all
          `>= 0`.
        - Prop. 3 (`delta = 4`): the analogous exact `F_4` (via the auxiliary
          `B_4 = sum a_j P_j^2`, `a_4 = 1/5`, `a_3 = (8b+5k+t-11)/20`, ...)
          and a mixed-square identity (9).  The bounds become positive after
          `t -> t+27`, `b -> b+11` or `k -> k+24`; the box minima are exact
          (4,764,410 quadratic minima).
        - Bridges: the atlas residual lists `(-1^5,+1,+1,-3,-3,-5)` (118) and
          `(-1^4,-2,-3,-4,-5)` (142); 324 signed multisets with labels `<= 4`;
          80 large-label lists.
        - FM-CHK67 (luna_max_venus, fresh code): ACCEPT Props. 2, 3 (all nine
          and twelve certificates, the 7,869-profile box, the finite
          reduction) and the coverage of every distance-3/4 list.
        - REPAIR (scope) of Prop. 1: the `-2 -> (+1, -1)` replacement applies
          only to NON-distinguished labels, and all counts refer to the
          labels other than the distinguished maximum.  Splitting a
          distinguished `-2` changes the distance: `(-2, -1, +1^7)` goes from
          distance 3 (EVEN 168, `F_3 = 84`) to 4.  With this scope,
          `EVEN = 2 F_delta` matched on 2,000 random lists (8 to 11 factors,
          distinguished labels up to 125), and clipping on 240/240.
      - **FM-MECH71 (astra_max_juno; `fm39/mech71_distance56_all_lists_repro.py`,
        rerun exactly in 48 s): FM3 holds for EVERY list at distance
        `delta <= 6`.  Also a tail theorem uniform in `delta`.**
        - Prop. 1: exact distance-5 and distance-6 polynomials through the
          core-multiset sums `C_(delta, sigma)` (formula (2) of the report).
          They use the shifted Krawtchouk `P_j` (`M = t - 2b`) and the
          multiplicities `M(q; sigma) = [U_0] U_1^q prod U_m`.
        - Certificates: matrix/scalar grid certificates (`d = 5`: 103,796 and
          4,520; `d = 6`: 1,405,291 and 112,523), and exact finite parts
          (`d = 6`: 22,593,136 quadratics, 67,779,911 points, minimum 0, i.e.
          support zeros only).  Also 1,101 small-list and 80 large-label
          character bridges.
        - Negative diagonal coefficients do occur, at positive words (136 at
          `d = 5`, 400 at `d = 6`), so the certificates are genuinely mixed.
        - Tail theorem (every `delta`):
          `max(t, b, k) >= (100 delta)^(300 delta^2)` implies `F_delta >= 0`.
          So every distance stratum is a finite exact computation.  The
          threshold is very coarse.
        - Independent check: FM-CHK69.
        - FM-CHK69 (luna_max_venus, fresh code).
          - Item 1 ACCEPT: formula (2)-(3), with 1,000 random compatible
            lists each at distances 5 and 6 (6-9 factors, distinguished
            labels up to 105) matching twice the formula value.
          - Item 2 REPAIR NEEDED (`d = 6` finite reduction).  The outer
            coefficient certificates pass.  The packet checks the displayed
            quadratic `Q_(y,z,w,q,eps)(i,j)` only at selected points.  It
            lacks a proof of `6 F_6 >= Q` for every admissible residual-pair
            allocation and all `i, j >= 0`, `i + j <= J`.
          - Item 3 REPAIR NEEDED (tail theorem).  The top-degree formula
            checks pass at `d = 5, 6`.  The general-`delta` majorants
            `E_j, O_j` (degree `<= delta - j - 1`, norms
            `<= (100 delta)^(100 delta^2)`, with the mixed signed-core and
            depletion terms) are not constructed in the packet.
          - Both repairs are delegated to the author as FM-MECH80.  Until
            then, distance 5 or 6 counts as proved with a pending repair.
          - Repairs (FM-MECH80, rerun: ALL CHECKS PASS).  Item 2: the exact
            gap identity (A2), `6 F_6 - Q(i,j) = 6[(q* - q_f) A_F +
            (q* - q_g) A_G + (q* - 1) A_H + q* A_L] >= 0`, holds for every
            residual-pair allocation.  Item 3: the general-`delta` majorant
            claim is withdrawn.  Nothing uses it: the tail statement now rests
            on FM-MECH77 (FM-CHK70).  The FM-MECH80 check is FM-CHK74.
      - FM-CHK65 (luna_max_eris, fresh code): ACCEPT all four items of
        FM-MECH63: the pair bounds (exact rational grids at 4,356 small, 660
        edge, 20,064 odd-label and 16,632 general edge points), the parity
        extraction, the resource table, the radial tails, and six exact
        words above the cutoff.
      - **FM-MECH69 (astra_max_vulcan; `fm39/mech69_many_labels_small_background_repro.py`,
        rerun exactly): any number of large labels on a labels-`<= 4`
        background.**  With `S = sum_i (n_i+1)^2` over the large labels
        (`n_i >= 5`, any signs and any number `k`), FM3 holds whenever
        `B >= max(384 S, 131072)`; also whenever `B >= 4096 S`.
        - New lemma: `|R_(n,-)|^2 <= 32 z` for every `n >= 5`.  Even `n` use
          `M - U_n(2q) >= (M/3)(1 - q^2)` and a second-derivative bound.  So
          two minus factors supply the needed `z` at bounded cost.  The
          stronger `|R/z| <= const` fails: `R_(4r,-)/z -> -(2r+1)`.
          Main-agent check: 50,000 random high-precision points, no
          violation (maximum ratio 0.19).
        - Exact residual after this theorem and the weighted theorem
          (`T >= 2^21`), with `Q = max(n_i + 1)` and
          `W = S + 16 J_3 + 25 J_4`: `W < 2^21 Q^2` and (`B < 384 S`, or
          `A < 256` and `B < 2^17`).  That is, few large labels in weighted
          count and a background small against them.  Ordinary counts are
          not bounded (e.g. `m` copies of 5 plus one label `m`).
        - FM-CHK68 (luna_max_neptune, fresh code): ACCEPT all five items.
          - The minus-factor lemma, via the sine pairing and the
            second-derivative bound (15,000+ exact profile checks).
          - The resource selection.
          - The radial constants (`81/524288`, `9/4096`, the `A < 256`
            branch); `k = 0` is assigned to FM-MECH48.
          - The exact residual (8).
          - Four exact words at `B = 131073` are positive.
      - FM-CHK66 (luna_max_pluto, fresh code).
        - ACCEPT FM-MECH61: the (BE) failure family, the T3/T4 margins, the
          cubic identity of Lemma 3, and the `omega < 64` screen (1,105,010
          pairs, minimum 67).
        - REPAIR (scope) FM-MECH66: the `e = 3` implication "(BE) at every
          long endpoint" must be restricted to the long-arc branch
          (`W < 0`).  At `(255,3,142,146)` the RF bracket is positive and
          (BE) fails, but `W > 0` there, so it is a short arc, covered
          separately.  Theorem 1 (`min(a,e) <= 3`) survives.  The checker
          also verified the `e = 3` RF formula on 48,104 pairs, the `e = 7`
          box, and a rational certificate for `theta/sin theta <= 11/8`.
      - FM-MECH70 (astra_max_ceres; `fm39/mech70_E_e456_repro.py`, rerun
        exactly in 27 s): (E), both signs, for every row with
        `min(a,e) <= 6` (PROVED), uniformly in the other exponent and both
        labels.
        - New cases: folded `e in {4,5,6}`.  The anchor caps, the lower bounds
          on `kappa` and the checkpoint constants are explicit for each `e`.
        - Lemma 2: the initial crossing satisfies (BE).  From the exact RF
          formula, a polynomial implication `G_e > 0 or V_e >= 0` is proved by
          Bernstein subdivision (37 terminal boxes at most, depth `<= 10`),
          with `S_8 <= exp`.
        - Lemma 3: the checkpoint gives `J` through the turning point, using
          the Sonin estimate and positive Taylor bounds.  Anchors beyond the
          cap are RF, since `K_e` has no root there.
        - Fixed box `sigma < 256`: 1,255,686 pairs, minimum about `4.8e13`.
          Also 3,380 RF bridges, 273 Catalan bridges and 43,383 large
          samples.
        - Remaining for (E) with `q >= 3`, both at `e >= 7`
          (`kappa >= 15`):
          - very central anchors `omega X < 6(sigma+1)`, `Y < omega/2`;
          - larger anchors `3 omega/5 < X <= sigma/2`,
            `X^2 < 4(e-1)(a+4)`.
      - FM-MECH67 (astra_max_minerva; `fm39/mech67_two_core_closed_forms_repro.py`,
        rerun exactly in 26 s): two core labels at every distance, partial.
        - Theorem (PROVED): both compatible signs for every two-core word
          with `n + m >= D = a + e + 2b`, at every distance and every
          `b >= 0`.
        - Prop. 1: a closed finite-binomial formula
          `A_k^(delta)(N,b) = sum_l V(delta,l) B(l,k)` for the diagonal
          coefficients (91 coefficients checked symbolically through
          `delta = 12`), and an explicit bilinear cross term (Prop. 3).
        - KILL (knob: coefficientwise positivity of the diagonal):
          - `A_(delta-2)^(delta)(2 delta - 2, 1) = -6(delta-4)/(delta(delta+1))
            < 0` for `delta >= 5`;
          - for every fixed `b`, the leading coefficient
            `[z^r] G(z)^b`, `G = 1 + 2z - 6z^2 + 2z^3 + z^4`, is negative at
            `r = 2b` (`b` odd) or `2b - 1` (`b` even), giving integer
            counterexamples at arbitrarily large distance.
          A Cauchy bound on the cross term also fails; a matrix repair works
          on one covered row only.
        - Remaining: two cores with `n + m < D`, `b >= 1`, distance `>= 6`.
          The adjacent layer `n + m = D - 2` is the first open one.
      - FM-MECH72 (astra_max_vulcan; `fm39/mech72_minus4_scaling_repro.py`,
        rerun exactly): a relative-error theorem in the residual, for a
        narrow sector.
        - Theorem: for one or two large labels `n_i >= 5` on a background of
          `h` factors `-4`, FM3 holds whenever `h >= max(S/8, 2^37)`,
          `S = sum (n_i+1)^2`.  The relative error against an explicit
          edge-plus-interior Gaussian model is `< 512/625`.  The model has a
          3D radial Gaussian at the edge and a 1D Gaussian at the interior
          maximum `y = sqrt(3/2)` of `|U_4(x) - U_4(y)|`.
        - Structural findings:
          - large distance does not force concentration at the corners
            `(+-2, +-2)`; edge and interior points both matter, and the model
            must keep the exact interior character values;
          - on FM-MECH53's no-margin family (`|F_d|/L_d -> 0`), the ratio to
            the corner model stays in `(9/10, 1)`.  So that obstruction was
            against the wrong model; a correctly scaled model may keep a
            margin.
        - Open: general labels-`<= 4` backgrounds in the residual (8).
      - FM-MECH73 (astra_max_ceres; `fm39/mech73_E_central_anchor_repro.py`,
        rerun exactly in 58 s): the very central anchors of (E), CLOSED,
        uniformly in `e`.
        - Lemma 2: the phase is fixed by reciprocity, using the centre
          coefficient for even `N` and `c_(-1) = (-1)^e c_1` for odd `N`.
          Each step turns by `theta` up to an error `arcsin(X_k/sigma)`.
        - Lemma 3: anchor-dependent crossing bounds on long arcs
          (`delta > 5/2` or `w > 13/5`, `53/20`, `21/10`, by range).  These
          give (BE) through a fixed certificate (degrees (11,10), 132
          coefficients, all `>= 0`), then Theorem C1 beyond.
        - With FM-MECH66, every anchor `X <= 3 omega/5` is covered.  Census:
          12,322 rows, 1,542,704 long pairs.
        - KILLs (knob: an unnecessary strengthening):
          - a uniform gap constant `delta >= 5/2` fails at
            `(103,10,58,62)` (where (BE) and (E) hold);
          - requiring `J(j,K) >= 0` at every larger anchor fails at
            `(160,9,120,125)`, but all 50 later endpoints there are short
            arcs satisfying (E).  So outer propagation needs a starting value
            only on the long-arc branch.
        - REMAINING for (E), `q >= 3`: larger anchors `3 omega/5 < X <= sigma/2`
          (forcing `omega/sigma < 5/6`), long arcs only, with
          `J >= 0` at `max(i_0, K)`.
      - FM-MECH75 (astra_max_vulcan; `fm39/mech75_h1_scaling_sectors_repro.py`):
        a relative-error theorem for the `H = 1` layer in two large sectors.
        - Model: `G = (2M/sqrt pi) int h(y) A(y)^(-3/2) exp(-M^2/(4A(y))) dmu(y)`,
          with `h(y) = (2-y)^e (2+y)^a (2+y^2)^b` and
          `A(y) = e/(2-y) + a/(2+y) + 4b/(2+y^2)`.  One angle is Gaussian and
          the other is integrated exactly; this keeps the quartic transition.
        - Theorem: with `M = p + 1 >= 2^10` and `l = floor(log_4 M)`, the
          bound `|F/G - 1| < 1/2` holds (so `F > 0`) in two sectors:
          - `b <= a, e <= 4b`, `b >= 2^16`, `3M^2 <= 16 l b`;
          - `r = min(a,e)`, `max(a,e) <= 4r`, `b <= r`, `r >= 2^22`,
            `M^2 <= 2 l r` (fixed `b`, even `b = 3`, allowed).
          These include words below the FM-MECH52 cutoff and the FM-MECH63
          `384 M^2` cutoff, with cutoff `O(p^2/log p)`.
        - No label-independent `D_0`; the `H = 1` layer stays open below these
          sectors.
      - FM-MECH74 (astra_max_ceres; `fm39/mech74_E_larger_anchor_repro.py`,
        rerun exactly).
        - Theorem 1 (PROVED): (E), both signs, at every endpoint for
          `3 omega/5 < X <= 7 omega/10` (`e >= 7`, `kappa >= 15`,
          `omega >= 64`), including `J` at `max(i_0, K)`.
          - A stronger crossing bound `kappa(y - x) > 2` on long arcs, from
            frozen binomial-energy coordinates.
          - (BE) at early long endpoints, by a (11,9) Bernstein certificate
            with 120 positive coefficients.
          - A checkpoint at `93 omega/100`, with margin 20691/160000.
        - Exterior direction: the boundary `c_(N+1) = 0` fixes the ratio
          `-u_k/v_k` for `X_k >= omega` through a monotone backward map.
        - Parameter-only test: if `4 s S <= pi^2 Delta^2` (with `s = K - j`)
          then the anchor has no long endpoint, so the short-arc theorem
          finishes it.  This covers the control `(160,9,120,125)`.
        - Remaining for (E), `q >= 3`: anchors `X > 7 omega/10` with
          `omega/sigma < 5/7`, outside the no-long-endpoint test, actual long
          endpoints only.  Sample positive residual: `(78,12,68,74)`.
      - **FM-MECH77 (astra_max_juno; `fm39/mech77_cubic_tail_repro.py`,
        rerun: ALL CHECKS PASS): a cubic tail threshold at every distance.**
        - Setting.  Distinguish a maximum label and split only the
          non-distinguished `-2` factors.  Then `t` = remaining label-1
          factors, `b` = `+2` factors, `k` = remaining labels `>= 3`,
          `sigma_3` = signed label-3 count, `N = t + 2b + k`, `s = delta + 1`.
        - Theorem: `N >= C_delta = 3200 s^3 + 760 s^2` implies FM3, for any
          labels, signs and background.  This replaces FM-MECH71's
          `(100 delta)^(300 delta^2)`.
        - Integer test (1): `N >= 380 s^2` and
          `N (3N - 1140 s^2)^2 >= 28800 sigma_3^2 s^3` implies FM3.  So
          `sigma_3 = 0` and `N >= 380 s^2` suffice.
        - Proof.  `EVEN = 2 F_delta`, with `F_delta` a coefficient of
          `(1 - z^2) E_eta prod f_(n,eps)` (Prop. 1).  The leading Gaussian
          form `L_d` is a positive sum of Hermite squares,
          `L_d = sum_j A(d,j) H_j^2` with `A(d,j) >= 0`.  Polynomial
          differentiation and multiplication bounds (Prop. 2) control every
          omitted logarithmic term: `|F_d - L_d| <= [(1 + rho^2) e^E - 1] L_d`
          with `E <= 228 (d+1)^2/N + 24 sqrt 2 |sigma_3| (d+1)^(3/2)/N^(3/2)`
          (Prop. 3).  Then `F_delta >= (133/864) L_delta > 0` (Prop. 4).
        - Prop. 5: the older uniform tests (weighted count, FM-MECH69,
          FM-MECH67) do NOT leave finitely many count profiles at any fixed
          `delta >= 7`.  Witness: `(+1)^m (-1)^m (+2)^3` with one label
          `2m + 6 - 2 delta` of sign `(-1)^m`, `m >= 2 delta + 10`.  The
          weighted test needs `delta >= 5,242,875` and FM-MECH69 needs
          `delta >= 13,822`.
        - Prop. 6: a finite representative family at each distance
          (labels clipped to `delta + 1`, `N < C_delta`).  At `delta = 7`,
          `C_7 = 1,687,040`.  The unfiltered box then has a 76-digit
          cardinality, far too large to enumerate.
        - Prop. 7 (kill, knob: fixed relative margin against `L_d` in a
          linear count region).  With `delta = 2r`, `a = e = cr` and the
          distinguished label `2(c-2) r`, `F/L -> 0`.  The family stays
          positive, so this does not exclude a linear FM3 threshold.
        - Independent check: FM-CHK70.  Open: lists with `N < C_delta` at
          `delta >= 7` (FM-MECH80 starts with `delta = 7`).
        - FM-CHK70 (luna_max_venus, fresh code): ACCEPT all four items.
          - The cubic threshold, and the fact that it does not use FM-MECH71's
            unconstructed general-`delta` majorants (FM-CHK69 item 3).  So
            the tail statement now rests on FM-MECH77 alone.
          - The integer criterion (1), including `sigma_3 = 0`.
          - Exact stress tests in the two stated families.
          - The finite representative description and the `delta = 7`
            bounds (the family itself is not enumerated).
      - FM-MECH76 (astra_max_minerva; `fm39/mech76_two_core_b1_band_repro.py`,
        rerun: PASS): two cores on a {1,2} background with `b >= 1`,
        `n + m < D`.
        - Combined rows `v_h(k) = [z^(k+h)] (1+z^2)^h P`.  The consumer is
          `R = T_b(j) - T_b(i) >= |W_b(j,i)|`.  The rows satisfy the coupled
          recurrence `(k+h+1) v_h(k+1) = Delta v_h(k) - (N+h-k+1) v_h(k-1)
          + 4h v_(h-1)(k)`.
        - `b = 1`: `W_1(j,i) = A_j B_i - B_j A_i` and
          `T_1(k) - T_1(k+1) = A_k B_(k+1) - B_k A_(k+1)`, so FM-MECH56's
          gluing identity survives.  A varying quadratic form keeps the
          exact square decrement
          `H_k - H_(k+1) = X (sigma c_(k-1) - Delta c_k)^2/(k+1)^2`.
        - Theorem 2 (PROVED): balanced `b = 1` (`a = e`), both signs, all
          labels, every distance.
        - Theorem 3 (PROVED): every `b >= 1`, both signs, on the band
          `4|Delta| <= sigma`, `sigma/4 + 2 <= n - m - 2b <= sigma/2`,
          `m >= 35 + 12b`.  The margin is explicit.
        - Failed (knob: monotonicity of the combined energy; unchanged Sonin
          rate; a two-dimensional determinant for `b >= 2`).
        - Open: the short-gap four-step inequality
          `T_1(j) - T_1(j+4) >= |A_j B_(j+4) - B_j A_(j+4)|` (`b = 1`,
          `m = 3`, `a != e`), then general short gaps (FM-MECH81).
      - FM-MECH81 (astra_max_minerva; `fm39/mech81_two_core_four_step_repro.py`,
        rerun: PASS): the four-step inequality
        `T_1(j) - T_1(j+4) >= |A_j B_(j+4) - B_j A_(j+4)|` for all `a, e >= 0`
        and `2 <= X = 2j - N <= N - 10`.
        - Lemma 2: with `p = c_j`, `q = c_(j+1)`, both signed slacks are
          binary quadratic forms in `(p, q)`.
        - Lemmas 3-5: uniform matrix regions, OL's central and outer ratio
          intervals, and an exhaustiveness split with no finite corner.  The
          two small folded levels use one polynomial with 82 positive
          coefficients.
        - Corollary 6: every two-core word with labels `n` and 3, exactly one
          `+2`, and any {1,2} background.  It uses FM-MECH71 for distance
          `<= 6` (two repairs pending) and FM-MECH67 for the boundary.
        - Matrix-only witnesses show that the row-direction information
          (Lemma 4) is needed.
        - Open: `b = 1`, `m >= 4` (next:
          `T_1(j) - T_1(j+5) >= |A_j B_(j+5) - B_j A_(j+5)|`), and `b >= 2`.
          FM-MECH84 asks for a mechanism uniform in `m`.  Check: FM-CHK73.
      - FM-CHK73 (luna_max_venus, fresh code) of FM-MECH76 and FM-MECH81.
        - ACCEPT: the combined-row consumer (both signs reproduce
          `2(R +- W_b)`; 80 direct Catalan cases, `b = 1, 2, 3`, labels up
          to 60).  Also FM-MECH76 Thm 2 (176,851 intervals through `L = 100`),
          FM-MECH81 Thm 1 (the 82-coefficient remainder rebuilt; the region
          split has no omitted corner), and Corollary 6 for exactly one
          `hat S_2`.
        - FM-MECH76 Thm 3: REPAIR (local, supplied by the checker).
          `omega >= 9 sigma/10` alone gives `25/(1296 sigma) > 1/(60 sigma)`.
          The hypothesis `4|Delta| <= sigma` gives
          `omega^2 = sigma^2 - Delta^2 >= 15 sigma^2/16`, which yields the
          cross-term bound `E_j/(60 sigma)` and `|T_b(i)| <= E_j/(184320 sigma)`.
          The stated margin follows, and the scope is unchanged.  All 826
          band rows through `N = 120` pass, plus 200 large rows.
      - FM-SEC140 (luna_max_jupiter; `fm39/sec140_atlas_update_repro.py`,
        phase 0 rerun exactly): the coverage atlas with FM-MECH63..70
        and ADV-2 added.  Phase 0 (sum `<= 30`, labels `<= 30`): 1,114,614
        lists, 123,678 uncovered.  Phase 1 (sum `<= 60`, labels `<= 40`,
        `h <= 3`): 78,037,215 lists, 10,059,259 uncovered.  No uncovered
        value is `<= 0`.  The smallest uncovered sum is 20, all at maximum
        label 5 and distance 5, e.g. `(-1^3, +1^2, +2, -3, -5, -5)` with
        value 76.  FM-MECH71 (distance 5 and 6) was not in this matcher.
      - **FM-MECH79 (astra_max_ceres; `fm39/mech79_E_last_anchors_repro.py`,
        rerun: PASS): the last (E) anchors.**  Both signs of (E) for
        folded rows with `a >= e+2`, `e >= 7`, `kappa >= 15`,
        `omega >= 64`, `7 omega/10 < X_j <= sigma/2`, at every endpoint.
        It includes `J(j,i) >= 0` at `i = max(i_0, K)`, which is the start
        of outer propagation.
        - Checkpoint estimate: `D_m <= r(x,z) D_j` with
          `r = [(1+z)/(1+x)] exp(-kappa (z^2 - x^2))`.  If `r < 1/100` and
          `w < 1/5`, the `J` margin is `>= 20691/160000`.
        - Crossing bound `kappa (y - x) > 14/5` on inner long endpoints,
          from the frozen binomial-energy coordinates.
        - Compact range `7/10 < x <= 49/50`: 28 bins `(k, Q, G)`.  For
          `kappa >= Q`, the bins use checkpoints at `G/10000` and a
          polynomial (BE) certificate with 6,272 Bernstein coefficients,
          all `> 1/25`.  For `kappa < Q`, a Jacobi-matrix root exclusion
          that is uniform in `a`: the normalized entries decrease with `A`.
        - Edge range `49/50 <= x < 1`: either `L = kappa h <= 11/(10 sqrt
          h)`, when the FM-MECH74 test `4 s S <= pi^2 Delta^2` excludes every
          long endpoint (`sS/Delta^2 < 12/5`); or a moving checkpoint
          `delta_0 = 9/(4 tau)`, `tau = h^(1/4)`, with `w < 1/10`.
        - With FM-MECH55..61, 66, 70, 73 and 74 this closes (E) for
          `q >= 3`, both signs, at every row and level.  So the two-label
          sector on the fundamental background is proved at every level.
          The region union is being checked independently (FM-CHK71).
        - FM-CHK71 (luna_max_eris, fresh code).
          - Item 1 ACCEPT: the proof steps of FM-MECH79.
          - Item 3 ACCEPT: the exact large-parameter spot checks.
          - Item 2 REPAIR (scope wording): state the anchor range
            `N/2 <= j`.  The consumer uses exactly this range (after folding,
            `j >= N/2`), and with it the union of FM-MECH55..61, 66, 70, 73,
            74 and 79 passes.  Adopted, so (E) holds at the consumer's
            strength.
      - FM-MECH78 (astra_max_vulcan; `fm39/mech78_h1_label10_repro.py` and
        `fm39/mech78_h1_constants_repro.py`, both rerun: PASS): `H = 1`
        constants.
        - The whole `H = 1` layer for labels `p <= 10`, every ratio and
          sign.  Bernstein-certified derivative bounds give the cutoffs
          `nu + b >= 4L + 2`: 166 for label 7, 214 for 9, 170 for 10+.
          Below them, 3,510,309 exact representatives remain after removing
          `b <= 2` and `d <= 6`.  All are positive; the least values are
          4,096 (7), 5,770 (8+), 4,850 (8-), 8,350 (9), 10,120 (10+) and
          10,348 (10-).
        - A positive-integral certificate valid for every ratio `a : e : b`
          (without the FM-MECH75 ratio restrictions).  Also a moderate-size
          `O(M^2/log M)` theorem on the quartic transition
          (`a = e = 2b`, e.g. `M = 17`, `b = 16`).
        - Kill (knob: a uniform relative-error margin of the FM-MECH75 model
          over all large distances).  The model has no such margin, so
          tightening its constants cannot finish `H = 1`.
        - Open, escaping every current criterion:
          `E[(x^2-y^2)^r (x^2+y^2-2)^3 (U_(2m^2)(x) - U_(2m^2)(y))] >= 0`,
          `r = m^2 + m - 3`, `m >= 23` (distance `m`).  This is passed to
          FM-MECH82.
        - FM-CHK72 (luna_max_neptune, fresh code): ACCEPT all three items:
          the closure through label 10 (cutoffs, reflection estimate,
          representatives), the ratio-free certificate, and the no-margin
          witnesses.
      - **FM-MECH80 (astra_max_juno; `fm39/mech80_distance7_all_lists_repro.py`,
        rerun: ALL CHECKS PASS in 284 s): FM3 for EVERY list at distance
        `delta <= 7`.**  It also repairs FM-CHK69 items 2 and 3 (recorded
        under FM-MECH71).
        - Prop. 3: an exact distance-7 polynomial with 39 nonempty kernels
          `C_(7,sigma)`.  It includes the mixed depletion terms
          `-(c_3 - 1) s_3 P_3 - (c_3 - 2) E_2(c_3, s_3)`: an unselected label
          3 is depleted when one or two other 3's are selected.  Two unsigned
          depletions first contribute at distance 8.
        - Prop. 4: 24 polynomial tail certificates give
          `t >= 249 or b >= 137 or k >= 333 => F_7 >= 0`, so `N >= 853`
          suffices.  The Schur determinant `D = K(R + z) - y^2` of the
          signed-count Gram matrix then removes 39,035 of the 39,615
          remaining rows, leaving 580 rows with `N <= 33`.
        - Prop. 5: in those rows the form is affine in the neutral-pair
          allocation, so its minimum is at one of five simplex vertices.
          8,157 small-corner allocations, 11,135,300 vertex quadratics and
          45,853,270 candidate evaluations give minimum `24 F_7 = 0`.
        - Knob: fixed `delta = 7`.  The next stratum `F_8` (`N < 2,394,360`
          outside the refined test) needs the two-depletion term.  Going one
          distance at a time does not reach the cone; what is needed is
          uniformity in `delta`.
        - FM-CHK74 (luna_max_neptune, fresh code).
          - ACCEPT items 1-3: the distance-6 gap identity, the distance-7
            polynomial (C1), the tail certificates and the Schur reduction.
          - Item 4 not yet established independently: the five-vertex
            argument checks, but the checker has not reproduced the 580-row
            set or its minima.  This is delegated as FM-CHK74b.  Until it
            returns, distance 7 counts as proved with that check pending.
          - FM-CHK74b (luna_max_neptune, fresh code): ACCEPT item 4.  The
            regenerated residual set has exactly 580 rows and matches the
            stated filter counts.  Direct enumeration of all admissible
            clipped allocations gives minimum `24 F_7 = 0`, with no negative
            row.  Distance 7 is closed and independently checked.
      - FM-MECH82 (astra_max_ceres; `fm39/mech82_h1_exterior_repro.py`,
        rerun: PASS): the `H = 1` consumer in combined rows.
        - `[U_p] g_(e,a,b) = T_b(k) - T_b(k+1)` with `k = (N+p)/2`.  The
          insertion identity is `T_(b+1)(k) - T_(b+1)(k+1) = T_b(k-1) -
          T_b(k+2) + W_b(k-1, k+2)`, so adding a `+2` needs the plus-sign,
          gap-three inequality on the previous background.  At `b = 0` this
          is Q2+.
        - Theorem (PROVED, uniform in `b`, `R` and `p`): `|a - e| <= 1`,
          `R = min(a,e) >= 1`, `b <= R`, `p = 2s + |a-e|`,
          `s^2 >= 4b(R+2)`.  Proof: every parity channel factors through
          Krawtchouk polynomials whose zeros are `<= 2 sqrt(h(R+2))`.
          This includes FM-MECH78's escaping family for every `m >= 4`.
        - Open: the gap-three step
          `T_b(k-1) - T_b(k+2) + W_b(k-1, k+2) >= 0` for `b >= 2` on
          general backgrounds, where the channel terms can be negative.
      - FM-MECH83 (astra_max_vulcan; `fm39/mech83_beta_saturation_repro.py`,
        rerun: PASS): two uniform regions.
        - Theorem 1: `F >= m(all)(1 - beta)`, with
          `beta = sum 1/(1 + max_(i in S) n_i)` over the negative contributing
          splits.  Proof: `U_a^2` contains every `U_(2j)`, `j <= a`, so
          `m(all) >= (1 + max_S n_i) m(S) m(S^c)`.  Consequence: for
          `L >= 4`, minimum label `>= 2^(L-2) - 3` gives FM3 at every
          distance and label ratio, with relative error
          `|V/2m(all) - 1| <= (2^(L-1) - L - 1)/(u+1)`.
        - Theorem 2 (exact saturation): with `h >= 2` non-distinguished
          labels above `delta` and the remaining background of weight
          `D <= min(delta + h - 2, 2 delta + 2)`, the value is
          `[q^delta] P(q)/(1-q)^(h-1) > 0`.  The proof uses the central
          moments of `P`, which are even and nonnegative.
        - Open: `h = 2`, `D > delta` (and `h <= 1`) in the residual
          `N < C_delta`.  An exact `r = 9` witness shows that any further
          argument must keep the helpful proper partitions.
        - FM-CHK75 (luna_max_venus, fresh code), together with FM-MECH82.
          - FM-MECH83 Thm 1 ACCEPT on the FM3 domain (even number of minus
            factors).  Odd-minus words vanish by `x <-> y`, and the `beta`
            expansion is not claimed for them.
          - FM-MECH83 Thm 2: REPAIR of the proof of (7), supplied by the
            checker.  The stated reason (a single selected large `U_n`
            vanishes by `eta`-degree) is false: background `U_2 U_2` has
            `<U_4, U_2^2> = 1`.  Replacement (split support): a contributing
            split `S` avoiding `p` has `w(S) <= 2 delta`, so it contains at
            most one label above `delta`, and none, since a large label
            would force `w(S) > 2 delta`.  The SU(2) coefficient formula
            then gives (7).  Coefficient identity, central moments and
            positivity pass.
          - FM-MECH82: REPAIR (convention).  In (4), `(L)_h` is the falling
            factorial `L (L-1) ... (L-h+1)`; with a rising factorial (4)
            fails at `R = 2, h = 2, L = 4, l = 0`.  With it, the coefficient
            bridge, the insertion identity, the exterior proof and the
            screens pass.  The negative channel at `(10,10,3,12)` (drops
            `[1575, 1400, 560, -546]` against a weighted total of 32214)
            refutes channelwise positivity only, not the consumer.
      - FM-SEC141 (luna_max_jupiter; `fm39/sec141_atlas_distance7_repro.py`,
        phase 0 rerun exactly): the atlas with FM-MECH71 (flagged) and
        FM-MECH77.  Phase 0: 100,743 uncovered lists; phase 1: 9,823,286.
        All uncovered values are positive.  The smallest uncovered sum is
        24, at distance 7, now covered by FM-MECH80.
        - Uncovered at distance 7 and 8: mostly the FM-MECH69 residual (8),
          i.e. few large labels on a small background (22,711 and 28,794 in
          phase 0).  The rest are two cores with `b >= 1` (995 and 1,133).
          There are none with `h = 1` and none with labels `<= 5`.
        - FM-MECH77 matches nothing inside the atlas ranges.
      - **FM-MECH85 (astra_max_juno; `fm39/mech85_quadratic_tail_repro.py`,
        rerun: ALL CHECKS PASS): a quadratic tail at every distance.**
        - Setting: `P(q) = E_eta prod_i (sum_(r <= n_i) q^r + eps_i q^(n_i/2)
          U_(n_i)(eta)) = sum_j p_j q^j` over the non-distinguished labels,
          so `EVEN = 2 F_delta` with `F_delta = p_delta - p_(delta-1)` (1).
        - Prop. 1: `p_j > 0` for every signed background.  The coefficient
          matrices are `I + eps D_n(c) J_n`, whose diagonal entries come
          from `Sym^n` of a rotation.
        - Theorems:
          - `N >= 432 (delta+1)^2` implies `F_delta > 0`, unconditionally,
            replacing FM-MECH77's cubic cutoff.
          - `k >= 40 delta` implies `F_delta >= (15/16) p_delta`.
          - One non-distinguished label above `delta` suffices, with any
            smaller background.  This removes FM-MECH83's support condition.
          - The general upper-band Gram/Schur formulas, and the complete
            depletion expansion.
        - Open (knob: an adjacent-coefficient comparison below the quadratic
          threshold): `p_delta - p_(delta-1) >= 0` on the explicit
          profiles (21)-(22), including the one-core layer `k = 0`
          (FM-MECH88 for `k >= 1`, FM-MECH89 for `k = 0`).  Check: FM-CHK76.
        - FM-CHK76 (luna_max_eris, fresh code): ACCEPT all four items (the
          quadratic tail with its constants, `k >= 40 delta`, Prop. 1 and the
          one-saturated-label statement, and exact tests at distances 8..14
          near the thresholds).
      - FM-MECH86 (astra_max_ceres; `fm39/mech86_h1_energy_repro.py`, rerun:
        PASS): an `H = 1` energy for every `b`.
        - The vector `V(k) = (v_0(k), ..., v_B(k))` satisfies a closed
          recurrence `t V(k+1) = A V(k) - u V(k-1)` with a constant matrix
          `A`.
        - The corrected energy `E_B(k) = S(|V(k)|^2 + |V(k-1)|^2) -
          2<V(k-1), A V(k)> + 4B T_(B-1)(k)` has the exact square decrement
          `E_B(k) - E_B(k+1) = (2k - N) |V(k-1) - V(k+1)|^2`.
        - The gap-three insertion step holds on an explicit region uniform
          in `b`, with unbounded labels, distances and numbers of `+2`.
        - Open: one comparison on the actual channel vectors.  Arbitrary
          vector states fail (exact witness), so the proof must use
          `v_(h+1)(k) = v_h(k-1) + v_h(k+1)`.  Delegated as FM-MECH89.
      - FM-MECH87 (astra_max_vulcan; `fm39/mech87_bounded_background_repro.py`
        (rerun: PASS, 353 s and 631 s) and
        `fm39/mech87_saturation_constants_repro.py` (rerun: PASS)): the
        dominant residual on bounded small backgrounds.
        - Reduction: with `T[a,j] = [z^a U_j(eta)] Q` and `c_i = T[2i,0]`,
          two large labels give `F = R + eps C` with
          `R = c_delta - c_(delta-n-1)`, `C = T[2 delta-n, n] - T[2 delta-n-2, n]`.
          Three large labels give an exact four-sign kernel.  For
          `delta >= D` positivity is analytic, and for `delta < D` the
          labels can be clipped to `delta + 1`.  So the finite search is
          `8 <= delta < D`, `max(5, 2 delta - D) <= m <= n <= delta + 1`,
          uniformly in the size of the large labels.
        - Certificates in exact 128-bit integers.  The absolute coefficient
          sum is `<= 4^60 < 2^127`, so nothing can overflow.  Scalar
          unimodality `c_0 <= ... <= c_(D/2)` holds for every background
          with `D <= 60`.  Two-label kernels: 1,307,953,746 (`D <= 60`).
          Three-label kernels: 4,321,391,262 (`D <= 48`).  No failure.
        - This covers all of FM-SEC141 phase 1, and its continuations with
          unbounded large labels.
        - Unbounded backgrounds: for `h >= 2` saturated labels,
          `delta >= ceil(15D/16)` suffices.  Also a count test
          `E = 16^(D-delta-1) R_0/lambda <= 1` via `P(1) >= lambda A_0` and
          absolute domination `sum |c_i| 16^-i <= A_0 R_0`, with an
          `h > 2` version.
        - Open: `h = 2`, `D > 60`, `delta < ceil(15D/16)` with `E > 1`, and
          `h <= 1` below saturation (FM-MECH90).
      - FM-MECH84 (astra_max_minerva; `fm39/mech84_two_core_short_arcs_repro.py`,
        rerun: PASS): two cores, `b = 1`, uniformly in the gap.
        - Index correction: the gap corresponds to `n = X + m` (the card's
          `X + m + 1` has odd degree).  The interior consumer includes the
          endpoints `i = N, N+1`.
        - Theorem (PROVED): for unbalanced rows (`a, e >= 1`) and
          `N/2 < j < i <= N+1`, if the polar sweep of `psi_j, ..., psi_i` is
          `< pi`, then `T_1(j) - T_1(i) >= |det(psi_j, psi_i)|`, with no bound
          on `i - j`.  The proof is a polygon argument from the one-label
          `b = 1` theorem applied to `(1 +- z) P`.
        - Also a normalized Sonin bound, an explicit root sector, and every
          gap with `min(a,e) = 1`.
        - Open: `b = 1` long arcs with `1 <= X_j < 2 + ((e-1)/(e+1)) omega`
          (a boxed sufficient criterion passes every exact test), and
          `b >= 2` (FM-MECH91).
      - FM-MECH89 (astra_max_ceres; `fm39/mech89_h1_criterion_repro.py`,
        rerun: PASS): `H = 1` from the actual channel relations.
        - Substituting the relations gives a block-tridiagonal formula for
          the insertion step, and an integer-only criterion (3), uniform in
          `B`, that does not assume positivity at `B - 1`.  It includes
          `n >= 5B` at `a = e` and `n >= 8B` with
          `|a - e| <= min(p, n/4)`.  These regions lie outside FM-MECH86's
          band, with unbounded backgrounds, labels and distances.
        - Kill (knob: treating consecutive odd-channel pairs as
          independent).  The unrestricted block Gram is false outside (3).
        - Open: the block-tridiagonal expression on the actual sequences
          outside (3), and the boundary range `p > N` (FM-MECH92).  Check:
          FM-CHK77.
      - FM-MECH90 (astra_max_vulcan; `fm39/mech90_h2_strata_repro.py`, rerun:
        PASS): regions with unbounded background weight.
        - Universal middle-distance identities, from 817 exact SOS identities
          and 234 central-coefficient identities.  These include an `h = 1`
          row: FM3 when the two largest labels have opposite signs and gap 0
          or 2, and the second largest is saturated, over any background.
        - A modulus tail certificate for every background.  Multiplying a
          complete `+n/-n` pair keeps cancellations that FM-MECH87's
          absolute bound discards.  The ellipse bound for `U_n(X)` is checked
          at 864 rational complex points.
        - Balanced backgrounds with `h = 2`: every `a, r >= 0` above an
          explicit boundary near `D/2`.  This is analytic for `H >= 44`; the
          990 profiles below that are checked exactly (116,589 inequalities,
          two implementations).
        - Kills (knobs: the modulus bound; individual-channel SOS).  At
          `r = 15` the channel sum is -756 while the full value is
          619,874,000,778.
        - Open: `F_r = [q^r] P_r(q)/(1-q) >= 0` with
          `P_r = E_eta[(G_3^2 - q^3 U_3^2)^r (G_4 + q^2 U_4)]`, and the
          left-of-centre `h = 2` residual for unbounded `D` (FM-MECH93).
          Check: FM-CHK78.
        - FM-CHK78 (luna_max_eris, fresh code): ACCEPT items 1-5 of FM-MECH87
          and FM-MECH90.  These are the two- and three-label kernels, the
          reduction to the finite search, independent exhaustive reruns at
          `D <= 40` (two labels) and `D <= 32` (three labels) with random
          samples to 60, the saturation band and the count test, and the
          FM-MECH90 strata.
      - **FM-SEC139 (main agent; `fm39/sec139_labels5_fast_repro.py`): the
        labels `<= 5` box is evaluated exactly.  All 227,336,512,420 entries
        are `>= 0` (3 zeros, least positive value 2).  With FM-SEC135, FM3
        holds for every word with labels `<= 5`.**
        - Evaluator.  Per background `(N, a)`, one mixed table
          `t[j][b] = E[s^(2a) d^(2(N-a)) P^j Z^b]`.  The core factors are
          in-place transforms: `hat S_4 = Z^2 + Z - 2P^2`, `Z - 1`,
          `h_4 = Z^2 + ZP - P^2 - 2P - 1`, `W_5 = hat S_5/s = Z^2 - ZP - P^2
          + 2P - 1`.  The label-3 factors enter by the exact combination
          `sum_j c_j t[j][b + l - j]`.  Everything is GMP integers, and
          every value is sign-checked.
        - Moments: an exact layer sweep up to degree 440.  The two long
          `N = 2` columns (to `b = 1013`) use the independent formula
          `sum_k c_k sum_i C(b,i) mu(k,i) mu(2N-k, b-i)`, and the two
          methods agree on 42 overlap points.
        - Checks in the run: 8 independent large-`b` moment values,
          342,468 shifted-moment bridges, and the direct Catalan values
          (300 polynomial words; the 181 in evaluated items all match).
          Per-item entry counts are checked against an independent
          enumeration, and the total against FM-SEC135's
          227,336,512,420.
        - Mirror symmetry.  `y -> -y` maps `s <-> d`, `P -> -P` and
          `h_4 <-> W_5`, which pairs background `a` with `N - a` for
          `1 <= a <= N-1`.  So only `a = 0` and `2a <= N` are evaluated, and
          the mirror half is counted.  The entry-count symmetry is checked
          exactly.  The values were checked directly on all 1,458 mirrored
          `(N, a, n)` groups with `N <= 12` (entries, zeros, least value).
        - Run: 384,374 work items `(N, a, n, v)`, at most 92 s each.  Each
          item is appended with fsync, and a rerun skips finished items.
          Kill-and-restart tests passed for both modes.  About 1.1 h at
          48 threads after the mirror restart.
        - Independent check: FM-CHK79.
        - FM-CHK79 (luna_max_pluto, fresh code): ACCEPT all five items.
          - The transform formulas and the mirror argument.
          - A fresh exact screen of 3,105 stratified box words: 44 at each
            `N = 2..71`, 21 large-factor records at `N = 20..40`, and four
            `N = 2` words with `b = 1002..1005`.  No negatives.
          - The three zeros `d^4 h_4`, `s^2 d^2 h_4` and `s d^2 S_5`, exact by
            two methods.
          - A code review found no indexing, dimension or overflow defect.
      - FM-MECH88 (astra_max_juno; `fm39/mech88_reflection_repro.py`, rerun:
        ALL CHECKS PASS): a neutral-pair reflection (2), an endpoint
        certificate (6) for arbitrary backgrounds, the families it gives at
        unbounded distance, and a terminal restriction (9) on the
        representatives.  A distance-30 example has a negative Gram pivot
        but positive `F`.
        - Open, `k = 1`: remove a non-distinguished core `n` of sign `eps`
          from background `B`; then `F_delta = R + eps C` with
          `R = b_delta - b_(delta-n-1)` and
          `C = T_B[2 delta-n, n] - T_B[2 delta-n-2, n]`.  Needed:
          `R >= |C|` on (9), starting at `n = delta - 1` (FM-MECH94).
      - FM-MECH93 (astra_max_vulcan; `fm39/mech93_saturated_sector_repro.py`,
        rerun: ALL CHECKS PASS): the saturated sector from FM-MECH85.
        - Every signed background has `c_j in Z_(>= 1)` for `0 <= j <= D`.
          The proof is an explicit positive-square construction, from
          `A_(n,eps)(c) = I + eps D_n(c) J_n > 0` with `D_n` the diagonal of
          `Sym^n` of a rotation.
        - Hence every `h >= 1` word is positive, with any background (no
          bound on `D`).  This includes `F_r >= r + 1` for FM-MECH90's
          family, whose single channel is -756 at `r = 15`.
        - What remains is the unsaturated case `h = 0`:
          `p_delta - p_(delta-1) >= 0` for `m <= delta <= floor((W-m)/2)`.
      - FM-MECH91 (astra_max_minerva; `fm39/mech91_two_core_outer_repro.py`,
        rerun: PASS): two cores, `b = 1`, long arcs.
        - Proved: uniform outer propagation of the `b = 1` criterion, and both
          FM3 signs on an explicit binomial tail sector (3).
        - Open: the inner long arcs `a > e >= 2`,
          `1 <= X_j < 2 + ((e-1)/(e+1)) omega`,
          `j + 5 <= i <= min(N+1, max(i_0, k_*))` with
          `k_* = ceil((N + omega)/2)` and
          `16 sigma^9 C(sigma, l+1) > C(sigma, j+1)`.  These need
          `2E >= tau H_j H_i` and
          `E^2 - tau H_j H_i E + C^2 (det Q)^2 H_j^2 H_i^2 >= 0`.
        - Missing supplier: a crossing/checkpoint estimate for the varying
          curve (6) (FM-MECH96).  The corrected energy alone gives no uniform
          endpoint bound at `b = 2`.
      - FM-MECH92 (astra_max_ceres; `fm39/mech92_h1_walk_expansion_repro.py`,
        rerun: PASS): an exact positive walk expansion for `H = 1`.
        - With `psi_j = (c_j, c_(j-1) + c_(j+1))` and the base-row chords
          `W(i,j) = det(psi_i, psi_j)`, the compatible channels are
          `v_h = K^h c`, where `K` is the shift sum.  Then
          `F_B(k) = [U_p] g_(e,a,B) = sum_(i<j) M_B(k;i,j) W(i,j)` with
          nonnegative integer weights `M_B` (weighted walks on odd-gap
          pairs), supported on `k - B <= i < j <= k + B + 1`.
        - Theorem (PROVED): with `t = min(a,e) >= 2`, the condition
          `(p - 2B)^2 >= 4(t-1)(N-t+2)` (3) gives both one-label families
          at every level and every `B`, including `p > N`.  Also the simpler
          linear-distance region `max(a,e) >= 15 min(a,e)`, `N >= 4d` (5).
        - Open: region (10), `(N - 2d)^2 < 4(t-1)(N-t+2)`.  (FM-MECH86/89
          regions: FM-CHK77 ACCEPT.)  The needed
          inequality is `(p+1) T_B(k) + |V(k)|^2 - |V(k+1)|^2 >=
          2B F_(B-1)(k)`.  Individual chords can be negative there.
        - Main-agent note: adjacent chords are Turan drops,
          `W(l, l+1) = D_l - D_(l+1)`, and the proved (E) bounds each chord by
          the drops between its endpoints.  A charging argument on the walk
          counts `M_B` may therefore close `H = 1`.  This is delegated as
          FM-MECH97.
        - Main-agent screens (`fm39/mech92_charging_screen.py`,
          `fm39/mech92_sweep_screen.py`).  The walk expansion matches the
          direct definition in 4,000 random cases (`a, e <= 40`, `B <= 8`),
          and every value is positive.
          - Kill (knob: using only the chord bounds from (E)).  The worst
            case over `|W(i,j)| <= D_i - D_j` is negative in 2,831 of 3,331
            right-half windows.
          - Window sweep `< pi` (so every chord `>= 0`) covers only 363 of
            4,364 right-half cases of region (10).  All chords are `>= 0` in
            635 of them.
          - So region (10) needs genuine cancellation between chords.
      - FM-SEC142 (luna_max_jupiter; `fm39/sec142_atlas_examples_check.py`
        reruns the exact example values and the table arithmetic; the census
        is a patched FM-SEC138 run).  The atlas with FM-MECH80 (flagged),
        FM-MECH81, 82 and FM-MECH83 Thms 1, 2.  Residuals: 75,559 (phase 0)
        and 9,583,360 (phase 1), all positive.  The smallest uncovered sum is
        26, at distance 8.
        - Not yet in this matcher: FM-SEC139 (labels `<= 5`), FM-MECH85
          (`h >= 1`, the quadratic tail, `k >= 40 delta`) and FM-MECH87.
          For example, the smallest witness `(14x(-1), +1, +2, -4, -5)` has
          labels `<= 5`, so it is covered.  FM-SEC143 adds them.
        - Distances 8-10: two cores with `b >= 1`, and the FM-MECH69
          small-background class.  At distance 8 the `h` profile is
          202,835 / 72,325 / 48 for `h = 0 / 1 / 2`; the `h >= 1` lists are
          now covered by FM-MECH85.
      - FM-CHK77 (luna_max_venus, fresh code) of FM-MECH86 and FM-MECH89:
        ACCEPT all items.  These are the closed recurrence, the corrected
        energy with its square decrement (and its sum-of-squares expansion),
        FM-MECH86 regions (5) and (7), and FM-MECH89 criterion (3) with both
        consequences.
      - FM-MECH95 (astra_max_vulcan; `fm39/mech95_double_endpoint_repro.py`,
        rerun: ALL CHECKS PASS): two sectors of the unsaturated `k >= 2`
        region.
        - Two non-distinguished labels equal to `delta`, any remaining
          background: `F = H_delta + (eps_1 + eps_2) A + eps_1 eps_2`, with
          `A = [z^delta] prod_B (1 + eps z^n)`.  The coefficient kernel gives
          `c_delta >= A^2/(delta+1)`, so `F >= (A + eps(delta+1))^2/(delta+1)
          >= 0`, or `F = H_delta - 1 >= delta` for opposite signs.  PROVED for
          every such word.
        - Quartet sector: with the four largest labels
          `p >= n_3 >= n_2 >= n_1` and the rest of weight `D`,
          `min(p, delta) - D + 1 >= 3 (floor(D/2) + 1)^2` implies FM3.  This
          replaces FM-MECH87's exponential loss by a quadratic one.
        - Open: the `k = 2` comparison (two cores `m <= n <= delta`, not both
          `delta`, on a {1,2} background), and the complement (12) for
          `k >= 3` (FM-MECH98).  Check: FM-CHK80.
      - Main-agent screen (`fm39/sec144_unsaturated_screen.cpp`).  The
          equivalent form is FM3(list) = `2 g_p`, with `g_p = [U_p(X)] E_y
          prod (U_(n_i)(X) + eps_i U_(n_i)(y))` over the background; this is
          checked against FM-SEC142's independent fusion DP on 200 random
          lists.  An exact GMP screen of the unsaturated region is running
          (see FM-SEC144).
      - FM-SEC144 (main agent; `fm39/sec144_unsaturated_screen.cpp`, run as
        `gsearch SAMPLES SEED THREADS 30 MODE`, and
        `fm39/sec144_h1_margin_scan.py B NMAX`): an exact screen of the
        unsaturated region.
        - 100,000 random exact values over five sampling modes, with
          `8 <= delta <= 30`: mixed; the edge `p = max`; many cores
          (`k <= 15`); H = 1 with many `+2`; two cores.  No negative value.
          The smallest `g_p / m_p(all)` per mode: 0.033, 0.031, 0.081,
          0.0055, 0.010.
        - H = 1 margin scan, exact and complete on the scanned ranges
          (background `(+1)^a (-1)^e (+2)^b`, distance `>= 8`, `p >= 6`):
          - With `a, e <= 26` and `b <= 12`, the minimum is always at
            `a = e = 26` and `p = 6`, falling smoothly by about 0.7 per extra
            `+2`, from 0.18 at `b = 0` to 0.0025 at `b = 12`.
          - At `a = e = 10` the minimum falls until `b` is about `a + 3`
            (0.023 at `b = 13`), then rises: 0.19 at `b = 25`, 4.6 at
            `b = 34`.
          - No negative value anywhere.
        - So the tight corner of H = 1 is balanced `a = e`, `b` near `a`,
          smallest `p`.  It lies outside the exterior theorem (FM-MECH82/92,
          which needs `p` large) and outside FM-MECH89 (`b <= a/5` at
          `a = e`).
        - Two-core scan (`fm39/sec144_two_core_margin_scan.py 8 16 9`,
          exact and complete): one extra core `+-m`, `3 <= m <= 9`, on
          `(+1)^a (-1)^e (+2)^b` with `a, e <= 16`, `b <= 8`, distance
          `>= 8`, unsaturated, `p >= max(m, 6)`.  No negative value.  For
          every `m` and sign, the minimum (0.015-0.018) is at balanced
          `a = e`, the largest `b` and the smallest admissible `p`.  So the
          hard corner is set by the background (balanced `+-1` with many
          `+2`), not by the core.
      - FM-MECH97 (astra_max_ceres; `fm39/mech97_h1_restricted_charging_repro.py`,
        rerun: PASS).
        - Kill (knob: charging from (E) alone, even with reciprocity and
          shift compatibility).  A uniform counterfamily exists, so the actual
          Krawtchouk recurrence is needed.
        - Restricted charging theorem (5)-(6) (PROVED): when every proper
          odd-gap chord in the window is `>= 0`, the outer chord is paid for
          by two overlapping subwindows, each sweeping `<= pi`.  This gives a
          new uniform H = 1 region for both parities, with no bound on `B`
          or the extra label.
        - Open: region (7), where some proper odd-gap chord in
          `[k-B, k+B+1]` is negative.  FM-MECH99 targets the tight corner of
          FM-SEC144 (`a = e`, `b ~ a`, small `p`).  Check: FM-CHK81.
        - FM-CHK81 (luna_max_neptune, fresh code).
          - FM-MECH88: ACCEPT (reflection (2), endpoint certificate (6),
            terminal reduction (9)).
          - FM-MECH97 restricted charging: REPAIR (statement).  State the
            identity for `B >= 1`, and handle `B = 0` separately by
            `F_0(k) = D_k - D_(k+1) = W(k, k+1) >= 0`.  The region then
            checks exactly.
          - The global-charging counterfamily: ACCEPT.  The charging bound is
            -9548 against the consumer value 5012.  An abstract row satisfies
            every tested (E) inequality but has `F_7 = -567`.
      - **FM-SEC145 (main agent): `+2` insertion monotonicity (conjecture,
        exact screens).**
        - Statement (M): for every signed background `B` and every
          `p >= max(labels of B, 3)`, `g_p(B + {+2}) >= g_p(B)`.
          Equivalently FM3(B, +2, p) >= FM3(B, p).  In the two-spin table
          `G(s,t) = sum_S eps_S m_t(S) m_s(S^c)` of `B` it reads
          `G(p-2, 0) + G(p+2, 0) + G(p, 2) >= 0`.  For H = 1, in FM-MECH82
          notation, it reads `F_B(k-1) + F_B(k+1) + W_B(k-1, k+2) >= 0`.
        - Consequence: removing the `+2` factors one at a time reduces every
          list to its `b = 0` version, with the same distinguished label.
          So (M) with Theorem OL proves the whole H = 1 layer, and (M) with
          (E) (FM-MECH79, FM-CHK71) proves the whole two-core layer, at every
          `b`.  It also reduces `k >= 2` to `b = 0` backgrounds.
        - Exact screens, no violation:
          - `fm39/sec145_h1_monotonicity_check.py 14 14`: H = 1,
            `a, e, b <= 14`, all `p >= 3`, 43,218 checks.
          - `fm39/sec145_z_insertion_monotonicity_screen.cpp` (GMP): random
            backgrounds with up to 10 cores (labels up to 40), up to 50 `+-1`
            and up to 12 `+2`, all `p >= max(labels, 3)`, 231,732 checks.
          - The smallest relative increase is 0, attained only on degenerate
            backgrounds (a single label equal to `p`).
        - Corner data (`fm39/sec144_corner_scan.py 6 40`): at `a = e`,
          `p = 6`, the minimum over `b` of `g_6/m_6(all)` is at `b = a + 3`
          and decays like `0.73^a` (`2.2e-6` at `a = 40`).  So arguments that
          perturb around the all-plus term cannot work there, while (M) is
          exact.
        - Delegated: FM-MECH100 (general proof) and FM-MECH101 (the H = 1
          case).
      - **FM-MECH102 (main agent; `fm39/mech102_pair_reduction_check.py`,
        PASS): the pair reduction.  FM3 for all lists follows from FM3 for
        PAIR-FREE lists, where each label value occurs with one sign only.
        Pair-free lists with at most one non-distinguished label `>= 3` are
        already proved.**
        - Identity:
          `(U_n(x) + U_n(y))(U_n(x) - U_n(y)) = U_n(x)^2 - U_n(y)^2
          = sum_(k=1..n) (U_(2k)(x) - U_(2k)(y))`.  So for every list with
          `+n` and `-n`, `FM3(rest, +n, -n) = sum_k FM3(rest, -2k)`.  Each
          summand has one factor fewer and the same minus parity; strong
          induction on the number of factors gives the reduction.  Checked
          symbolically for `n <= 12` and on 2,000 random lists.  Examples:
          `(+1)(-1) = (-2)`, `(-2)(+2) = (-2) + (-4)`.
        - Pair-free H = 1 (no non-distinguished core): the background is
          `(+-1)^a (+-2)^b` with one sign per label.
          - With `-2` present or no 2's, there is no `+2`, so Theorem OL
            applies (the `-2 = ds` split is allowed).
          - With `+2` present, the 1's share one sign, so `min(a,e) = 0` and
            FM-MECH52 Thm 1 applies.
        - Pair-free two cores (one non-distinguished core `m >= 3`):
          - without `+2`, (E) on the fundamental background (FM-MECH79,
            FM-CHK71);
          - with `+2`, `min(a,e) = 0` and the FM-MECH64 strip theorem
            (every distance, every `b`, all labels) applies.  FM-MECH64 has
            been rerun but not yet independently checked; that check is
            FM-CHK82.
        - So the whole cone reduces to pair-free lists with at least TWO
          non-distinguished labels `>= 3`, within the unsaturated region and
          outside the other proved strata.  Non-pair-free H = 1 and two-core
          lists decompose into such lists (reducing `(-2)(+2)` creates `-4`
          cores).  For example, FM-SEC144's tight corner `a = e`, `b ~ a`
          is not pair-free.
      - FM-SEC146 (main agent; `fm39/sec146_pair_free_screen.cpp`,
        `fm39/sec146_pair_free_adversarial.py 240`): margins of the pair-free
        residual (`k >= 2` non-distinguished cores, unsaturated, distance
        `>= 8`).
        - Random exact screens: 50,000 lists.  No negative value; the
          smallest `g_p / m_p(all)` is 0.94 with few 1's and 0.999 with many
          single-sign 1's.
        - Adversarial hill-climbing (32 restarts, 4 min each, weight
          `<= 90`): no negative value; the smallest ratio is 0.805.  The
          extremal shapes are nearly all-minus with many small cores and
          `p` near the maximum, e.g. background
          `(-8, -7^4, -6^2, -5^3, -4^2, -3^3, -2, -1)` with `p = 9`.
        - Atlas residual shapes are far from tight.  For
          `(22 x (-1), -8, -8; +22)`, `g_p / m_p(all) = 39,782`.
        - Contrast: the non-pair-free H = 1 corner decays like `0.73^a`
          (FM-SEC144).  So pair-free lists looked robustly positive, with
          `g_p >= 0.8 m_p(all)` in every case seen.
        - CORRECTED (00:10): there is NO uniform relative margin on
          pair-free lists.  FM-MECH103 (Juno) and FM-MECH109 (Vulcan) give
          explicit pair-free families, including full all-minus lists at
          bounded weighted count, with `g_p / m_p(all) -> 0` (below
          `1e-80`).  The adversarial search above did not reach them.
          Arguments on the pair-free residual must therefore be exact.
        - Larger adversarial run (`sec146_pair_free_adversarial.py 420
          150`, weight `<= 150`): no negative value, and the smallest ratio
          is 0.758.  Structured all-minus families `(-3..-m)^r`
          (`sec146_pair_free_families.py 7`, `m = 5, 6, 7`, `r <= 7`): ratios
          0.95-1.26, settling as `r` grows, with the minimum at `delta = 8`.
          The pair-free margin shows no sign of tending to 0.
      - FM-MECH94 (astra_max_juno; `fm39/mech94_two_core_band_repro.py`,
        rerun: ALL CHECKS PASS): two cores, the band `n = delta - 1`
        (`p + n = D - 2`).  PROVED for every {1,2} background and both
        signs, so coverage extends to `p + n >= D - 2`.  The next band
        (`p + n = D - 4`) needs a coupled estimate.  After FM-MECH102 the
        two-core interior is no longer on the critical path.
      - FM-MECH100 (luna_max_mars): `+2` insertion monotonicity (M) for
        every signed background when `p >= ` the total background weight
        (the weight-dominant H = 1 and two-core subregions).  The interior
        `max(3, max n_i) <= p < sum n_i` is open.  No standalone verifier
        was supplied, so this is recorded as a claim only.
      - FM-SEC147 (main agent; `fm39/sec147_mono_*.py`): insertion
        monotonicity on pair-free backgrounds (conjectures, exact screens).
        - (M+): adding a plus label to a pair-free background (even `n`
          singly, odd `n` in pairs) never decreases `g_p`, for
          `p >= max(labels, 3)`.  No violation: `+4`, `+6`, `+8` (16,538,
          17,107, 18,775 checks), `(+1,+1)`, `(+3,+3)`, `(+5,+5)` (7,840,
          14,379, 16,516 checks), plus FM-SEC145's `+2`.
        - (M-): inside all-minus backgrounds (single odd plus labels
          allowed), adding minus labels in pairs never decreases `g_p`.  No
          violation: an existing class `(-n,-n)` (114,212 checks), a new
          class `(-n,-n)` (118,597), two new distinct labels `(-n,-m)`
          (95,755).
        - Kill (knob: minus insertion into mixed backgrounds).  Adding
          `(-3,-3)` to the all-plus background `(1^12, 2^4, 4, 6, 6)` lowers
          `g_6`.  So plus labels must be stripped before minus labels.
        - Also two distinct odd plus labels `(+n,+m)`: 78,933 checks, no
          violation (`sec147_mono_plus_nm.py`).  C++ GMP stress test
          (`sec147_monotone_screen.cpp`, both lemmas, up to 10 cores, labels
          up to 30, many 1's): 389,874 checks, no violation.
      - **FM-MECH106 (main agent; `fm39/mech106_monotone_reduction_check.py`,
        PASS): FM3 for the whole cone follows from two insertion-monotonicity
        lemmas on pair-free backgrounds.**
        - (M+): for a pair-free background `B` and `p >= max(labels, 3)`,
          adding plus labels (one even label, or two odd labels) keeping
          pair-freeness never decreases `g_p`.
        - (M-): for an all-minus background with at most one odd plus label,
          adding two minus labels never decreases `g_p`.
        - Descent: from a pair-free list (FM-MECH102 reduces every list to
          these), strip even plus labels one at a time and odd plus labels
          in pairs by (M+).  Then strip minus labels in pairs by (M-), keeping
          `p`.  `p <= 2` is the proved {1,2} sector.
        - The base lists have at most three factors and are `>= 0`: `(p)`
          gives 0, `(+-n, +-p)` gives `2[n = p]`, and `(+o, -n, -p)` gives
          `2 m(o,n,p)`.  The script checks 14,470 base values (labels
          `<= 40`) and replays 1,180 descent steps on random pair-free lists.
        - Status: (M+) and (M-) are conjectures (FM-SEC147: about 900,000
          exact checks, no violation).  `+2` insertion (FM-SEC145, FM-MECH100,
          FM-MECH101) is the `n = 2` case of (M+).  This route needs neither
          H = 1 nor two-core results.
      - FM-CHK83 (luna_max_neptune, fresh code): FM-MECH64's strip theorem
        is ACCEPTED, with one wording refinement.  `s` preserves the full
        antisymmetric cone, while `Z` preserves its equal-parity subcone.
        Items 2-4 (consumer, exact tests, scope) pass.  This confirms the
        pair-free two-core `+2` component used in FM-MECH102.
      - Returned at 23:10; partial, verifiers not yet rerun by the main agent
        (recorded as claims).  The H = 1 and two-core items are now off the
        critical path (FM-MECH102/106).
        - FM-MECH99 (astra_max_ceres): balanced H = 1 for
          `L >= s(s+1) + 1` (`a = e = b = L`, `p = 2s`).
        - FM-MECH96 (astra_max_minerva): two cores, `b = 1`, both signs for
          `3 omega/5 < X_j <= sigma/2` and for all `omega < 64`.
        - FM-MECH98 (astra_max_vulcan): the complete `k = 2` first interior
          strip; uniform cutoffs `delta >= 32(t+2)` and `delta >= 64(t+3)`;
          an arbitrary-background resource theorem.  Also: `m_p(all)`
          admits no fixed positive relative margin even at `k = 2`.
          Main-agent check: the family uses `(+1)^r (-1)^r (+2)^r` with
          cores `(-6q, +7q; -13q)`, which is NOT pair-free (`+1` and `-1`
          both occur).  So the pair reduction removes it, and it does not
          conflict with FM-SEC146's pair-free minimum of 0.758.
        - FM-MECH101 (luna_max_mercury): the H = 1 insertion step in an
          exterior region.  Kill (knob: termwise channel positivity).
        - FM-MECH104 (luna_max_mars): an exact kernel sum with a ballot-walk
          formula for all-minus backgrounds.  Kill (knob: termwise kernel
          positivity); the aggregate stays positive.
      - **CORRECTION (main agent, 23:45) to FM-SEC147 and FM-MECH106: route A
        is killed as stated.**
        - Screening error.  When an inserted or removed pair has odd total
          weight, the smaller background has `g_p = 0` by parity at the same
          `p`, so the comparison is trivially true.  It is also useless for a
          descent, since it assumes what it should prove.  This made vacuous:
          the `(-n,-m)` screen (`sec147_mono_minus_nm.py`) whenever `n + m`
          was odd; mode 1 of `sec147_monotone_screen.cpp` likewise; and the
          "two largest" test, which counted odd-sum removals.  The other
          FM-SEC147 screens (`+even`, `(+odd, +odd)`, `(-n,-n)`) have even
          weight and are not affected.
        - Kills (exact, even-weight).
          - (M-) as stated: inserting `(-2,-2)` into the all-minus
            `(-5, -3^6, -1^5)` lowers `g_6` from 713,552 to 243,743.  Other
            witnesses: `(-8,-2)`, `(-2,-4)`.
          - (M+) for odd pairs: inserting `(+1,+1)` into
            `(-9^5, 4^6, 12^3)` lowers `g_13`; also `(+3,+3)`, `(+1,+3)`,
            `(+5,+1)`.  This holds even when the pair is the two largest odd
            plus labels.
          - Single even MINUS insertions: `-4` into
            `(1^5, 2^5, 5^3, 12^4)` lowers `g_12`.
          - The majority-by-weight removal rule on minus and odd-plus
            backgrounds: removing `(-12, -2)` from
            `(-12, -9^6, -5, -3^6, -2, -1^5)` raises `g_12`.
        - Survivors, with even weight and adversarial search, no violation:
          - (M+even), one even PLUS label (`fm39/sec147_mplus_adv.py ... even`);
          - odd plus pairs into backgrounds with at most one minus label
            (T2).
        - Conjecture (D): for every pair-free `B` and admissible `p`, SOME
          removal (one even label, or two labels of even total weight) does
          not increase `g_p`.  It holds on 56,542 random backgrounds; an
          adversarial test is pending.  Without a rule (D) does not give a
          proof.
        - FM-MECH106's logic (descent plus base lists) stays valid, but its
          hypotheses (M+) and (M-) are false as stated.  Route B (pair-free
          `k >= 2` directly) and possible rule-based descents remain.
      - FM-CHK80 (luna_max_venus, fresh code): ACCEPT FM-MECH84 (`b = 1`
        short arcs, root and minimum-label sectors), FM-MECH91 (`b = 1` tail
        sector and propagation from a valid seed), FM-MECH92 (H = 1 exterior
        sector including `p > N`) and FM-MECH95 (double-endpoint sector and
        quartets with (5)).
      - FM-CHK82 (luna_max_eris, fresh code): FM-MECH102 ACCEPT (identity,
        induction, and the `k <= 1` case analysis), with FM-MECH64 cited in
        its equal-parity wording (FM-CHK83).  No negative value.
      - FM-CHK84 (luna_max_neptune, fresh code): REJECT (M+) and (M-) as
        stated.  The failures are exact, pair-free and consumer-compatible,
        and agree with the 23:45 correction; none is a negative FM3 word.
      - FM-MECH103 (astra_max_juno; `fm39/mech103_pair_free_no_margin_repro.py`,
        rerun: ALL CHECKS PASS).
        - Kill (knob: a uniform relative margin `g_p >= c m_p(all)` on
          pair-free lists).  An explicit all-minus family with weighted
          count `T < 903/400` has ratio below `1e-80` (at `q = 50`); an exact
          finite example already beats 0.8.
        - Also killed: scalar log-concavity of the coefficients, and
          positivity of the quadratic form on arbitrary recurrence states.
        - Open: the all-minus `k = 2` inequality `F_d >= 0` (8) on the actual
          row `c_j = [z^j](1+z)^b (1-z)^(a+b)` (FM-MECH114).
      - FM-MECH109 (astra_max_vulcan; `fm39/mech109_family4_repro.py`, rerun:
        ALL CHECKS PASS).  FM3 for every word
        `(-p)^2 (-6)(-8) prod_(j=1..10) (-2j)^(2 l r_j)` uniformly in `l`,
        over its whole unsaturated range of `p`.  Also: no uniform positive
        margin on pair-free lists, including strictly unsaturated lists
        with weighted count below 3.
      - FM-MECH105 (luna_max_pluto; `fm39/mech105_sign_budget_repro.py`,
        rerun: PASS): a sign-counted maximum-label budget.  Sort
        `3 <= n_1 <= ... <= n_L`; let `c_j` count the subsets `A` of `[j]`
        with `j in A` and an odd number of minus labels, and
        `C = sum_(j < L) c_j/(n_j + 1)`.  Then FM-MECH83's bound gives
        `F >= M([L])(1 - C)`, so `C <= 1` proves FM3.  This covers every
        list of at most four factors, all-plus lists, and the pure-core
        pair-free region `C <= 1`.  Open: `C > 1`, via a toggle inequality
        on lower lists (FM-MECH116).  Check: FM-CHK85.
      - FM-MECH107 (astra_max_ceres; `fm39/mech107_mplus_dominant_minus_repro.py`,
        rerun: PASS): (M+) on the dominant-minus-core class.  The
        background is `s^a Z^b (U_n(x) - U_n(y)) prod_i (U_(m_i)(x) +
        U_(m_i)(y))` with `n >= M` and `a >= eta` (4).  The proof is uniform
        in the number and size of the plus cores and in `a, b`.  The
        one-background-core `+2` step follows.  Open: the two-core ordering
        `B = (+1)^a (+2)^b (-m, +n)`, `m < n` (FM-MECH117).
      - FM-MECH110 (luna_max_mercury; `fm39/mech110_mplus_h1_repro.py`, rerun:
        PASS): (M+) with `+2` for every pair-free H = 1 background
        `(eps 1)^a (+2)^b`, and the stated character-valued core extensions.
        The FM-MECH82/92 criterion needs `min(a,e) >= 2`, so this is a new
        direct proof.  Open: arbitrary signed cores (FM-MECH118).  Check:
        FM-CHK85.
      - FM-SEC143 (luna_max_jupiter): the atlas with FM-SEC139, FM-MECH85
        (`h >= 1`, tails) and FM-MECH87.  Phase 1 (sum `<= 60`, `h <= 3`)
        is entirely covered.  Phase 0 (sum `<= 30`) leaves 1,203 lists, all
        positive.  The filtered frontier starts at `delta = 8` with `k = 3`
        non-distinguished cores, and `k = 4` appears among the first 20.
        The verifier is a patched matcher run, not rerun by the main agent.
      - FM-MECH113 (luna_max_mars): exact identity for even `n`,
        `g_p(B + {+n}) - g_p(B) = G(p,n) + sum_(j=1..n/2) (g_(p-2j)(B) +
        g_(p+2j)(B))`.  So (M+even) reduces to bounding the cross entry
        `G(p,n)` by shifted coefficients.  Proved for backgrounds of length
        `<= 1`, and screened exactly on signed (also non-pair-free)
        backgrounds of length `<= 4` with labels `<= 4`.  Recorded as a
        claim; the general estimate is open.
      - FM-MECH118 (luna_max_mercury; `fm39/mech118_mplus_genuine_repro.py`,
        rerun: PASS).  For an even insertion `2q <= p`,
        `g_p(B (+2q)) - g_p(B) = sum_(j=1..q) (G_B(p-2j, 0) + G_B(p+2j, 0))
        + G_B(p, 2q)`.  So (M+even) holds whenever `B` is a genuine
        `SU(2) x SU(2)` character (all `G_B >= 0`), e.g. all-plus
        backgrounds and the pair-free H = 1 slice.  Exact two-core screen:
        no failure.  Open: the same two-core ordering as FM-MECH117.
      - FM-MECH119 (luna_max_mars).  Kill (knob: the POINTWISE shifted
        coefficient kernel `K^Delta_ij(c) = K_ij(c) - K_(i-1,j-1)(c)` of
        FM-MECH85).  For the pair-free all-minus `(-4)^20` with `p = 40`,
        `delta = 20`, the value is `F_20 = p_20 - p_19 = 113,693,884,711`
        (direct fusion gives `2 F_20`).  At `c = 99/100`, `K(c)` has 21
        positive pivots, but `K^Delta(c)` has a negative pivot at index 18.
        The integrated (over `c`) version is not excluded (FM-MECH120).
      - FM-MECH115 (astra_max_vulcan; `fm39/mech115_pair_free_k34_repro.py`,
        all six cases rerun: PASS): every PAIR-FREE list with three or four
        non-distinguished cores at distance 8, 9 or 10, with any number of
        single-sign 1's and 2's.
        - `+2` sector: a Newton expansion
          `F_delta(a,b) = sum C_ij C(a,i) C(b,j)`.  At `A = 2 delta + 4` the
          verifier proves `Delta_a^i Delta_b^j F(A,0) >= 0` and
          `Delta_b^j F(a,A) >= 0`, plus the finite box `a, b < A`.
        - `-2` sector: `f_(2,-) = (1+z^2+eta z)(1+z^2-eta z)` reduces the
          background to fundamentals.  An exact identity
          `(2 delta)! F = sum_r D_r(t) K_floor(r/2) K_ceil(r/2)` with monic
          Krawtchouk `K_j` and a coefficientwise dominance (4) beyond a
          threshold, plus finite parts.
        - Finite values (k, delta): (3,8) 154,036; (3,9) 308,474; (3,10)
          586,806; (4,8) 545,216; (4,9) 1,175,300; (4,10) 2,433,900.  All
          `>= 0`, with independent two-spin controls.
        - Knob: fixed `delta`, since the boxes grow with `delta`.  This does
          not reach the cone by itself; FM-MECH121 asks for uniformity in
          `delta`.  Open next: `k = 5` at `delta = 8` and `(k, delta) = (3, 11)`.
      - FM-MECH116 (luna_max_pluto; `fm39/mech116_prefix_screen.py`): the
        toggle route.  Pick a minus factor `i`, let `R` be the rest, and let
        `T_s` be the sum of the toggle terms over subsets `A` of `R` of size
        `s`, so that FM3 is `sum_s T_s`.
        - Kills: every layer `T_s >= 0` (at the `delta = 8` frontier
          `T_2 = T_3 = -1`); every single toggle determinant `>= 0` (a `-1`
          from `A = (-3,-4)`); and `(+1,+1)` insertion (`g_13` drops by
          9,603,276,419).
        - Surviving conjecture (P): every cumulative prefix
          `sum_(s <= t) T_s` is `>= 0`.  It passes Pluto's screens, but has
          not been tested adversarially (FM-MECH122).
      - FM-MECH120 (luna_max_mars).  The c-integrated kernel
        `bar K_ij = sum_h 2/(i+j-2h+2) C(i+j-2h, i-h) v_(i+j-2h,h)`
        (measure `2c dc`) has `bar K_dd = p_d`, so the consumer is the corner
        `(bar K - shift)_dd = p_d - p_(d-1)`.  Its shifted version is PSD on
        `(-4)^20` and on the tight Jupiter and Juno backgrounds.
        - Kill (knob: positivity on every feature direction).  On
          `(-3)^5 (-5)^6` (`delta = 20`) the shifted kernel has a negative
          pivot exactly at index 20, while `F_20 = 4,598,165 > 0`.  So a PSD
          certificate of the shifted kernel is stronger than FM3 at the
          consumer index itself.
      - **FM-SEC148 (luna_max_venus) and FM-SEC151 (main agent): a monotone
        descent rule.**
        - FM-SEC148, exhaustive: every pair-free background of total weight
          `W <= 32` (835,218 backgrounds with `>= 2` factors; 8,839,302
          admissible `(B,p)`; 67,964,688 comparisons).  Conjecture (D)
          holds throughout.  Rule 1 has zero failures there: remove two
          copies of a most frequent label (ties: smallest label); if all
          multiplicities are 1, the largest same-parity pair; else the
          smallest even label.  It does not depend on `p`.  Three
          alternative rules fail, with smallest witness
          `(-1)^15 (-2)^2` at `(W,p) = (19,3)`.
        - FM-SEC151 (`fm39/sec151_rule1_counterexample.py`): Rule 1 FAILS at
          `W = 97`.  For `B = (-1^2, 2^2, 6^2, 7^2, 8^2, 9^2, 10^2, 11)` and
          `p = 11`, all multiplicities tie, so the rule removes `(-1,-1)`,
          and `g_11` rises from 4,589,803,663 to 6,291,848,748.  Twenty-two
          other removals are monotone there, including every single even
          plus label.
        - Rule 1' = remove an even plus label first (smallest), else Rule 1.
          Adversarial hill-climbing finds no violation: 32 runs of 300 s at
          `W <= 100` (`sec151_rule1p_adversarial.py 300 100
          evenplus_first`), and 32 runs of 420 s at `W <= 140` biased to
          backgrounds without even plus labels and with near-tied
          multiplicities (`sec151_rule1p_partb_adversarial.py 420 140
          partb`).  The only zeros are trivial backgrounds.
        - LEMMA R (conjecture): the Rule 1' removal never increases `g_p`.
          It splits into (R-a) = (M+even) and (R-b), the duplicate-first
          removal on backgrounds without even plus labels.  Since every
          Rule 1' removal has even weight and lowers the number of labels,
          Lemma R with FM-MECH102 would give FM3 for EVERY list: descend to
          a background with at most one label, where `g_p >= 0` trivially.
      - Returns at 00:55 (verifiers rerun by the main agent: all PASS).
        - FM-MECH112 (astra_max_minerva; `fm39/mech112_quartet_descent_repro.py`):
          an exact quartet translation law, and a descent that preserves the
          pair-free consumer (odd plus factors allowed in the fixed
          background), with explicit infinite positive all-minus families.
          Kill (knob: complete fusion-channel positivity).
        - FM-MECH114 (astra_max_juno; `fm39/mech114_no1_subcone_repro.py`):
          all-minus `k = 2` with no label-1 factor (`a = 0`) and at least one
          odd non-distinguished core, by parity fusion and a binomial
          certificate.  Also `delta = 8`, `k = 3` with maximum 6 or 7, by
          exhaustive enumeration.  Kills: continuous relaxations of the
          Krawtchouk row.
        - FM-MECH117 (astra_max_ceres; `fm39/mech117_two_core_plus2_repro.py`):
          the reversed two-core `+2` comparison for every `a, b >= 0` and
          `3 <= m < n < p`, by exchanging the distinguished label with the
          minus core before the dominant-minus argument; regions (8) for
          even insertions with many plus cores.
        - FM-CHK86 (luna_max_mercury): ACCEPT FM-MECH107 on its gap-restricted
          cone (proof inequalities, parity and gap conditions, the finite
          box).  The unrestricted cone action is false (exact
          counterexample), so the statement keeps the gap restriction.
        - FM-SEC150 (luna_max_mars): the FM-MECH115 method at `k = 3` through
          `delta = 13`, with no failure.  Profile and box sizes grow
          polynomially; the `-2` thresholds grow in steps.
        - FM-SEC149 (luna_max_jupiter): the pair-free residual has 330 words
          through sum 30 and 80,212 through sum 40, all positive (680 roots
          after the pair-reduction tree in phase 0).  First pair-free
          example: `(-1)^6 (-5)^2 (-6)^2`, value 4146.
      - FM-MECH121 (astra_max_vulcan; `fm39/mech121_linear_thresholds_repro.py`,
        rerun: PASS): delta-uniform linear resource thresholds for PAIR-FREE
        lists, at every `delta >= 2` and independently of the core labels
        (`N = a + 2b + k`):
        - `k = 3`: `N >= 22 delta - 2` implies FM3;
        - `k = 4`: `N >= 36 delta - 1` implies FM3;
        - a general-`k` version.
        Also an explicit same-gap Krawtchouk identity, and mixed Newton
        positivity at shifts `4 delta + 2` (`k = 3`) and `5 delta + 3`
        (`k = 4`).  Open: the shift `2 delta + 4` and the boxes below the
        thresholds, which grow with `delta`.  A fixed positive fraction of
        `P_delta` is unavailable (exact endpoint family).
      - FM-CHK85 (luna_max_eris, fresh code).
          - ACCEPT FM-MECH110's H = 1 (M+) proof.
          - ACCEPT FM-MECH105: the budget theorem, the all-plus case and the
            finite screens.
          - REPAIR NEEDED (a side claim only): FM-MECH110's bounded
            core-census count omitted repeated core labels and
            double-counted the `a = 0` sign choice.  The proof does not
            depend on that count.
      - **FM-SEC152 (luna_max_venus; `fm39/sec152_rule1p_exhaustive.py`,
        rerun by the main agent: PASS): Rule 1' (Lemma R) is exhaustively
        verified for every pair-free background with total weight
        `W <= 40`.**  That is 6,011,040 backgrounds with `>= 2` factors and
        81,808,160 admissible `(B,p)` cases, with zero failures (W <= 32:
        8.8M cases; W <= 36: 27.8M).  Arithmetic is exact signed 128-bit
        with the bound `prod 2(n_i+1) <= 4^W`.
        - Smallest relative slack per branch: even-plus step 0 (trivially,
          at `B = (+2,-3)`, `p = 3`); duplicate step 7/4 (at
          `B = (-1^2, -2^2, -3)`, `p = 3`, removing `(-1,-1)`: 11 versus 4);
          fallback 0 (trivially, at `B = (-2,-3)`).  No child value is
          negative.
        - The printed verifier had a duplicated struct block (a
          transcription error).  The fm39 copy deletes the duplicate and
          reproduces the published counts exactly.
      - **KILL (01:15): LEMMA R (Rule 1') is FALSE as stated.**  Exact
        counterexamples (FM-MECH124 astra_max_minerva, FM-MECH125
        astra_max_juno; checked independently by the main agent with its own
        two-spin recursion, `fm39/sec155_lemmaR_counterexamples.py`):
        - `B = (-1)(-2)^2 (-3)^2 (-5)^2 ... (-15)^2` (all-minus, `W = 131`),
          `p = 15`.  The rule removes `(-2)^2`; `g_15` goes from
          703,572,078,442 to 1,087,917,863,463.
        - `B = (-1)^2 (+3)^2 (+5)^2 ... (+15)^2` (`W = 128`), `p = 16`.  The
          rule removes `(-1,-1)`; `g_16` goes from 214,784,436,431 to
          296,201,364,755.  Reflecting the odd signs gives a witness with a
          uniquely most frequent label, so it is not a tie artifact.
        - Knob: the prescribed frequency-based selection.  In both witnesses
          another removal is monotone (for example `(-15)^2` in the first),
          so the existence form (D) is not refuted.  FM-SEC152's exhaustive
          range `W <= 40` and the adversarial search to `W = 140` missed
          these structured families (every label doubled, a run of odd
          labels).
        - FM-MECH122 (luna_max_pluto): REJECT conjecture (P) (cumulative
          prefixes).  It fails for the pure-core `(4^24, -)` word at every
          choice of the minus factor.
        - FM-MECH126 (luna_max_mars): the distinct-label step of Rule 1'
          (steps 3-4) is proved for lengths 2, 3 and 4 (claim, not rerun).
        - FM-CHK87 (luna_max_mercury): ACCEPT FM-MECH114 (the parity-fusion
          subcone and the `delta = 8`, `k = 3` frontier for maximum 6 or 7)
          and FM-MECH117 (the reversed two-core `+2` comparison).
        - Status: the cone is NOT proved.  It reduces (FM-MECH102) to
          pair-free lists with `k >= 2`.  The open descent statement is now
          (D): some even-weight removal is monotone for each `(B, p)`.  (D)
          is exhaustive for `W <= 32` (FM-SEC148), but no selection rule is
          known.
      - **FM-SEC156 (main agent): Rule W, a weight-based descent rule (new
        conjecture LEMMA W).**
        - Diagnosis (`fm39/sec156_D_all_removals.cpp`).  In every Lemma R
          witness (W = 97, 128, 131, 362) nearly every even-weight removal is
          monotone for ALL `p` at once.  The rule's choice was one of the few
          failing ones, always the duplicate of the SMALLEST label.  (D) also
          holds on the structured families up to `W = 885`.
        - Rule W.  Take the label class of largest total weight
          (multiplicity times label; ties go to the larger label).  If it has
          two or more copies, remove two.  Otherwise pair it with the largest
          other label of the same parity, or remove it alone if it is even.
          Every removal has even weight.  REPAIR (FM-CHK89): if the heaviest
          class is a single odd label with no odd partner, remove the
          smallest even label (step W4).  Example: `B = (-3,-2)`.  The
          exhaustive C++ sweep already used this fallback.
        - Evidence, no failure:
          - exhaustive over every pair-free background with `W <= 40`
            (81,808,160 `(B,p)` cases; `fm39/sec156_ruleW_exhaustive.py`,
            Venus's sweep with Rule W);
          - every Lemma R witness and the structured families `F1`, `F2`
            (odd runs with multiplicities 2, 3) up to `W = 2030`, every `p`
            (`fm39/sec156_ruleW_single.cpp`);
          - adversarial hill-climbing, 32 x 360 s at `W <= 120`, seeded
            with flat multiplicities (`fm39/sec156_ruleW_adversarial.py`).
            The smallest relative increase is 1.0.
        - LEMMA W (conjecture): the Rule W removal never increases `g_p`.
          With FM-MECH102 it would prove FM3 for every list.  Each step is
          even-weight and lowers the number of labels, down to at most one
          label.
      - Returns 01:20 (claims; verifiers not rerun):
        - FM-MECH127 (astra_max_vulcan): uniform comparison thresholds for
          both small-label removals; new three- and four-core sectors with
          large labels.
        - FM-MECH123 (astra_max_ceres): (M+even) on two-minus backgrounds
          (`+2` with `a >= b`; every even insertion on bare two-minus).
        - FM-SEC153 (luna_max_jupiter): the Rule 1' descent on all 80,212
          (sum `<= 40`) and 365,281 (sum `<= 44`) pair-free residual roots
          never fails, consistent with Rule 1' failing only above `W = 97`.
      - FM-MECH132 (luna_max_mars): exact reductions of Lemma W's
        single-copy branches (claim, not rerun).  With
        `G_C(s,t) = sum_S eps_S m_s(C - S) m_t(S)` and fusion channels
        `I(a,b)`:
        - W3 (remove an even singleton `eps n`) is equivalent to
          `sum_(s in I(p,n), s != p) G_C(s,0) >= |G_C(p,n)|`, for all labels
          of `C` odd and every class weight of `C` at most `n`;
        - W2 (remove `(eps m, eta n)`) is an explicit signed comparison in
          `G_C`.
        Both are proved for all-plus `C`.  The knob is mixed-sign
        cancellation.
      - FM-CHK89 (luna_max_neptune, fresh code): an adversarial search for
        Lemma W violations up to `W = 200` (164 structured seeds and 32
        hill-climb restarts, exact) finds none.  The smallest normalized
        slack is 1, at `B = (-1)^32`, `R = (-1,-1)`, `p = 32` (support
        edge: `g = 1` versus 0).  It also found the undefined branch
        `B = (-3,-2)`, now repaired by step W4.
      - **KILL (02:40): LEMMA W is FALSE**, by three independent exact
        counterexamples.  The main agent checked the first with
        `fm39/sec156_ruleW_single.cpp`.
        - FM-MECH128 (astra_max_minerva): the all-minus
          `B = (-1)^10 (-8)^4 prod_(v = 3, 5, ..., 31) (-v)^floor(31/v)`
          (`W = 427`).  Class 8 is uniquely heaviest (weight 32 against at
          most 31), so the rule removes `(-8)^2`.  At `p = 31`, `g_31` rises
          from 35,976,660,898,477,897,403,704,068,501,225,333,918,277,825 to
          36,075,697,823,243,267,121,589,443,884,412,616,296,612,401.  This is
          the only failing `p` of 199.
        - FM-MECH129 (astra_max_juno): a mixed background with a unique
          heaviest class, `A = (+1)^35 (+3)^11 (+5)^7 ...`, failing for both
          signs of the selected odd label.
        - FM-MECH131 (astra_max_ceres): a W1 decrease with a uniquely
          heaviest class of eight `+9` (weight 72).
        - (D) still holds on the `W = 427` witness: every pair containing
          the maximum label `-31` is monotone for all `p`
          (`fm39/sec156_D_all_removals.cpp`).
      - KILLS (main agent) of other selection and averaging candidates:
        - Rule M (by the maximum label: two copies if repeated, else with
          the largest same-parity label, else alone if even, else the
          smallest even label; `fm39/sec159_ruleM_exhaustive.py`).
          Exhaustive through `W = 40`, it fails 234 times even in the
          unproved region (`p >= 6`, `delta >= 8`, unsaturated); first at
          `B = (-1)^22 (-2)^2`, `p = 6`, removing `(-2,-2)`.  Rule W handles
          that case, and Rule M handles `W = 427`.
        - Averaged (D) with uniform weights over all even-weight removals:
          negative at some `p` for mixed backgrounds near `W = 140-160`
          (adversarial).
        - Label-weighted (`w_R = prod n_i`), dimension-weighted
          (`w_R = prod (n_i + 1)`) and dimension-squared-weighted averages were
          nonnegative on every earlier witness
          (`fm39/sec159_weighted_averages.cpp`).  Adversarial search then
          KILLS them too (`sec159_weighted_averages_adversarial.py`):
          - dimension-weighted: negative at 3 `p` on
            `(-1^2, 3, 4^5, 5^7, 6^7, 7^2, 8, 16^2)`;
          - label-weighted: negative at 5 `p` on
            `(-17^3, -5, -3^6, 1^2, 2^2, 4^7, 6^8, 12, 16^2)`.
          Every rule and averaging scheme tried has failed.  The bare
          existence statement (D) has survived every test, but it gives no
          proof without a selection principle.
      - FM-CHK88 (luna_max_eris, fresh code): ACCEPT FM-MECH112 and
        FM-MECH121.
      - FM-MECH133 (luna_max_pluto): Rule W holds for all-plus backgrounds
        and for parity-aligned signs `eps_i = (-1)^(n_i)`.  Main-agent note:
        the reflection `y -> -y` flips the signs of odd labels and leaves
        `g_p` invariant, so the parity-aligned stratum is the mirror of the
        all-plus one.
      - **FM-SEC160 (main agent): conjecture (S), a two-child averaging
        inequality on the residual.**
        - Residual region: pair-free `(B, p)` with `p >= max(labels)`,
          `p >= 6`, `delta = (W - p)/2 >= 8`, unsaturated (`max(B) <= delta`)
          and `k >= 2` labels `>= 3` in `B`.  Its complement is proved:
          labels `<= 5` (FM-SEC139, FM-CHK79); distance `<= 7`
          (FM-MECH68/71/80, checked); saturated `h >= 1` (FM-MECH85,
          FM-CHK76); pair-free `k <= 1` (FM-MECH102, FM-CHK82/83); and
          `p <= 2` (the {1,2} sector).
        - (S): `2 g_p(B) >= g_p(B - R_W) + g_p(B - R_M)` on the residual.
          Here `R_W` is the Rule W removal (FM-SEC156 with step W4) and
          `R_M` the Rule M removal (the maximum label: two copies if
          repeated, else with the largest other same-parity label, else
          alone if even, else the smallest even label).
        - Why it suffices: induct on the number of labels.  The children
          have fewer labels and are either residual (induction) or proved,
          so both are `>= 0` and `g_p(B)` is at least their average.  With
          FM-MECH102 this would give FM3 for every list.
        - Evidence:
          - exhaustive over every pair-free background with `W <= 40`: zero
            failures in 22,358,566 residual `(B, p)` cases
            (`fm39/sec160_S_exhaustive.py`; the 132 raw failures all lie
            outside the residual, e.g. `(-1)^17 (-2)` at `p = 3`);
          - the four known single-rule counterexamples (`W = 34`, 427, 540,
            869), every `p`: the two rules never fail together and (S)
            holds (`fm39/sec160_WM_single.cpp`).
          - Adversarial test in the residual (`fm39/sec160_S_adversarial.py
            600 260`): (S) is KILLED.  On `B = (-16, -13^6, -9^3, -7^3,
            -3^3, -2, -1^11)` (`W = 164`), (S) fails at 9 residual `p`; on
            `B = (-10, 7^12, 9^8)` (`W = 166`) at 2.  In both, Rule W alone is
            monotone at every `p` and (D) holds; Rule M's failure there is too
            large for the average.
          - LEMMA WM (conjecture): on the residual, at least one of
            `R_W`, `R_M` is monotone, i.e.
            `max(g_p(B) - g_p(B - R_W), g_p(B) - g_p(B - R_M)) >= 0`.  With
            FM-MECH102 and the proved strata, this gives FM3 for every list
            by induction.  Evidence: never both failing exhaustively through
            `W <= 40`, and on all seven known single-rule counterexamples.
            Adversarial search finds the minimum of
            `max(Delta_W, Delta_M) / g_p(B)` to be 0.648 (`W <= 260`) and
            0.665 (`W <= 700`), with no violation
            (`sec160_W_or_M_adversarial.py`, `..._large.py`).  The margin is
            robust; a proof needs a criterion for which rule applies.
        - Context: Rule W alone fails only at large `W` (Minerva 427, Juno
          540, Ceres 869), and is clean through `W = 48` exhaustively
          (FM-SEC157, Venus) and on 1,368,356 residual roots (FM-SEC158,
          Jupiter).  Rule M alone fails 142 times in the residual through
          `W = 40` (first `B = (-1^4, -2, -3^8, -4)`, `p = 6`), but handles all
          three large Rule W failures.
      - Returns 03:10 (claims, not rerun): FM-MECH130 (astra_max_vulcan):
        Rule W holds for every background with at most four factors and in
        unsaturated three- and four-core regions.  FM-CHK90
        (luna_max_mercury): ACCEPT FM-MECH123; FM-MECH127 needs a repair
        (its uniform radial estimates are not derived).  FM-MECH134
        (luna_max_mars): W3 for `|C| <= 2` and an all-minus parity invariant.
      - FM-CHK91 (luna_max_neptune, fresh `cpp_int` code): no Rule W violation
        on 40 pair-free backgrounds with `208 <= W <= 600` (5,780 `p`
        checks, including Juno's `(-1)^2 (+3)(+5)...(+47)`).  This is a small
        sample: the known Lemma W counterexamples (W = 427, 540, 869) were
        not in it, so it does not affect the kill.
      - FM-MECH140 (luna_max_mars): Rule M is monotone on all-plus and
        parity-aligned backgrounds, at every `p >= max`, by nonnegative
        Clebsch-Gordan coefficients (claim, not rerun).  On mixed signs, Rule M
        alone fails already at the max-class weight ratio
        `rho = mult(max) max / max_v (v mult(v)) = 1/6`, where Rule W works.
      - Main-agent sweep (`fm39/sec160_WM_exhaustive48.py`, running):
        Lemma WM has zero failures over ALL pair-free `(B, p)` through
        `W = 40` (81,808,160 cases), not only on the residual.
      - **KILL (04:25): LEMMA WM is FALSE** (FM-MECH137, astra_max_juno;
        checked independently by the main agent with the pruned single-`p`
        evaluator `fm39/sec162_gp_single.cpp`).
        - `B = (-11)^5 (-54)(-2) A`, where `A` holds each positive odd label
          `n` in `1..53` except 11 with multiplicity `floor(54/n)`.  Then
          `W = 1283`, `p = 55`, 135 factors, `delta = 614`, `k = 80`, so
          this is a residual case.
        - Rule W uniquely picks `(-11,-11)` (weight 55), and Rule M picks
          `(-54,-2)`.  Both increase `g_55`: from 8.4993e119 to 1.1488e120 and
          to 1.3591e120.
        - The removal `(+1,+1)` is monotone (5.33e118), so (D) holds.
        - Pattern: every fixed candidate set (Rule 1', Rule W, Rule M, W or
          M) has failed at some larger weight.  The working removal moves
          between the smallest labels, the largest labels and the heaviest
          class.  (D) itself has never failed.  The lemma's exhaustive range
          (zero failures over all pair-free `(B,p)` through `W = 48`,
          `fm39/sec160_WM_exhaustive48.py`, finished 04:35) is far below this
          witness.  Finite exhaustive ranges therefore say little about these
          selection lemmas.
      - FM-SEC164 (luna_max_mars; `fm39/sec164_D_exhaustive.py`, rerun by
        the main agent at `W <= 36`: PASS).  (D) is exhaustively verified on
        the residual through `W = 40`: every residual `(B, p)` (22,358,566
        pairs, 222,096,442 removal comparisons) has a monotone even-weight
        removal.  The margin is large: the minimum over `(B,p)` of
        `max_R (g_p(B) - g_p(B-R)) / g_p(B)` is `22449/25937` (about 0.866), at
        `B = (-1^8, -2^4, -3^2)`, `p = 6`, `R = (-1,-3)`.  Every residual
        parent value is positive.
      - FM-CHK92 (luna_max_eris, fresh code).  ACCEPT FM-MECH130's
        Rule W result for at most four factors, on the range
        `p >= max(3, max |B_i|)` (state this restriction), and its three-
        and four-core results and channel obstruction.  REPAIR (wording):
        the `(-1,-4)` witness lies outside the consumer `p`-range.
      - FM-MECH144 (astra_max_juno; `fm39/mech144_many_ones_D_repro.py`,
        rerun: ALL CHECKS PASS): (D) and FM3 on a uniform many-ones region.
        - Reflect odd signs so the fundamentals are `+1`, and write
          `B = (+1)^a C` with no label 1 in `C`.  Let `t` be the number of
          minus factors in `C`, `l = ceil(t/2)`, `nu = 2l + 3`, and
          `K = 1 + sum_(+n in C) 5n(n+2) + sum_(-n in C) 3(n-1)(n+3)` plus
          `5p(p+2)` (`t` even) or `3(p-1)(p+3)` (`t` odd).
        - If `Q(a) = 15 a^2 - K nu (2a + nu) - 15 >= 0`, then
          `g_p((+1)^a C) >= g_p((+1)^(a-2) C)` and `g_p > 0`.  This holds
          for arbitrary core labels, signs and factor count, and gives both
          (D) and FM3 on that region.
        - Also exact certificates that tell removable core pairs from
          non-removable ones.  Open: the low-fundamental-count residual.
          The W = 1283 witness lies outside the uniform bound; the verifier
          prints its core threshold `a >= 276920`.
      - FM-MECH143 (astra_max_ceres; `fm39/mech143_two_odd_fusion_repro.py`,
        rerun: PASS).  Exact positive identities.
        - TWO-ODD FUSION (proved).  Let `C` be any signed list of even
          labels and `a, b` odd.  Then `Phi(C, eps a, eta b) = sum_(j in
          CG(a,b)) Phi(C, (eps eta) j)`, where `+0` counts `2 Phi(C)` and
          `-0` counts 0.  Reason: `C` is even in `x` and in `y` separately,
          so the two cross terms integrate to zero.  Every child has fewer
          factors.  With FM-MECH102, a minimal counterexample to FM3 (by
          factor count) is pair-free and has either no odd labels or at
          least four.
        - Level-2 sector (proved): `Phi(-a,-b,-c,-d,+2^t,+e_1..+e_s) >= 0`
          for odd `a, b`, even `c, d`, even `e_i >= 4`, every `t`, whenever
          `d >= a + b + sum e_i`.  The proof uses an ordered antisymmetric
          cone `A_ij = U_i(x)U_j(y) - U_j(x)U_i(y)` with the gap condition
          kept (not the cone action rejected in FM-CHK86).
        - Obstruction, knob "one block fusion with a nonnegative
          remainder".  `Q_I = 2^(1-m) sum_(J subset I, |J| even)
          Phi(Lambda^J)`, so a block fusion lowers the value exactly when a
          sign flip inside the block does.  At `(-1,+2,+3,+4,-5,+6,+7,+8)`
          (residual with `p = 8`, `W = 28`, `delta = 10`), `Phi = 956` is the
          minimum over all even sign patterns, so all 247 block fusions
          increase it (960..1310).  (D) holds there: `(-5,+7)` gives
          `g_8 = 18 <= 478`.  This does not refute FM3 or (D).
        - A context-free polynomial identity with nonnegative coefficients
          into fewer factors or lower weight is impossible for an all-minus
          product (lowest order in `h` at `x = z + h, y = z`).  Averaged
          identities remain available.
        - Consequence for minimal counterexamples: `Phi(Lambda) +
          Phi(Lambda^(ij)) = 2 Q_(ij) >= 0` for every pair by induction.  So
          a negative `Phi(Lambda)` forces every double flip to be at least
          `|Phi(Lambda)|`.
        - Open (first unresolved formula): four odd labels on an even
          background, `Phi = Q_1 + Q_2 + Q_3 - 2H`; this needs
          `Q_1 + Q_2 + Q_3 >= 2H`, which is equivalent to FM3 there.
      - FM-MECH142 (astra_max_minerva; `fm39/mech142_separated_quartet_repro.py`,
        rerun: PASS; coefficient stress test `fm39/mech142_hc_stress.py` by
        the main agent: 608 boundary cases with `D < 60`, `M <= 15`, minimum
        `h_c = 4`).  Separated-quartet descent (proved).
        - Let `C` be any signed word of weight `D`, and append a pair-free
          quartet `n_i = q + t_i`, `0 = t_1 <= ... <= t_4 = M`, with
          `q >= D + M + 1`.  Then `Phi(C, quartet)` equals a sum of
          `Phi(C, alpha_P a, beta_P b)` over the three pair partitions
          `P = I|J`, `a in CG(n_I)`, `b in CG(n_J)`, `a, b > 0`, `a + b <= D`,
          plus `sum_c h_c Phi(C, sigma c) + 2 h_0 Phi(C)`, with integer
          `h_c = m_c + z_c - N_c >= 0`.  The 3|1 splits vanish by degree.
          Every child has at least two fewer factors; children may have
          mixed signs and pairs (FM-MECH102 repairs those).
        - This proves every compatible pair-free 7- or 8-factor word that
          satisfies the separation (with the six-factor sector), and in
          general reduces separated words to shorter ones.  Knob: quartet
          separation.
        - Remark (main agent): the same degree argument for two labels gives
          `Phi(C, e1 n1, e2 n2) = sum_(c in CG(n1,n2)) Phi(C, (e1 e2) c)`
          whenever `n1 + n2 > D`.  With `n1 = p` and `n2 = max(B)` this is
          exactly `max(B) > delta`, the proved `h >= 1` stratum.  So the
          residual (`max(B) <= delta`) is exactly where the two-label
          version fails.
        - Kills (knob: unrestricted quartet translation, including a shift
          by 2): lowering the quartet increases the value, at
          `(-1)^12 (-6)^3 (-8)` (11,273,148 -> 12,382,358) and at
          `(-1)^26 (-2)(-3)^2 (-4)(-7)^4` (shift 2).  FM3 is not affected.
        - No uniform projected identity: an exact functional `L` is
          nonnegative on every compatible word with at most three factors
          after projection to `a + b <= 18`, but `L = -62,922` on
          `V_7- V_8- V_9- V_10-`.  Knob: using only the degree bound of `C`.
          Identities that use the actual coefficients of `C` remain open.
        - Proved family: every all-minus word of eight consecutive labels
          `(-r)...(-(r+7))` (closed quintic polynomials with positive
          coefficients for `r >= 14`, exact values for `r <= 13`).
      - FM-MECH145 (astra_max_vulcan; `fm39/mech145_delta_uniform_k34_repro.py`,
        rerun: PASS).  Route B, delta-uniform regions of the pair-free
        `k = 3, 4` layers below the linear cutoffs of FM-MECH121 (fundamental
        signs reflected to `+1`; `a` ones, `b` twos; all core signs).
        - `+2` background: `a >= 2 delta + 4`, `delta >= 2048`, every core
          label `>= ceil(delta/2)` (`k = 3, 4`; exact Newton certificate
          `X <= 1`, minimum core label of order `sqrt(delta log delta)` in
          the integer form); `b >= 2 delta + 4` with core labels `>= 6`
          (`k = 3`) or `>= 7` (`k = 4`), any `a`.
        - `-2` background: `a + 2b >= 2 delta`, `delta >= 8192`, core labels
          `>= ceil(delta/2)`; and for `k = 3`, `b >= 1`, `a <= 4b`,
          `delta >= 65536`, core labels `>= ceil(delta/2)` with no further
          count restriction (the balanced part of the finite boxes).
        - Exact complement stated in the report (finite for each `delta`,
          `O(delta^(k+2))` profiles).  Knobs: a fixed positive relative
          margin cannot cover all core profiles (three `+3` cores at
          `a = 2 delta + 4` have `F_delta/P_delta -> 0`); in the `-2` sector
          the top mixed coefficient `2 Cat(delta-1) - Cat(delta)` is negative
          (`-572` at `delta = 8`).
      - FM-SEC163 (luna_max_venus; `fm39/sec163_selector_census.py`,
        `fm39/sec163_witness_masks.py`, `fm39/sec163_min_selector_sets.py`;
        witness masks rerun by the main agent and identical; the W <= 40
        census not rerun, its (D) count agrees with FM-SEC164).
        - Over all 22,358,566 residual `(B,p)` through `W = 40`, every case
          has 1 to 19 monotone removals.  Of 12 fixed selectors only Rule W
          covers the whole range; Rule M misses 142 cases, the largest
          allowed pair 40, the top two compatible factors 261,978.
        - With the eight large witnesses (W = 128, 131, 34, 427, 540 and its
          reflection, 869, 1283) there are eight minimum two-selector
          portfolios, e.g. {Rule W, duplicate most frequent}.  These are
          finite-screen portfolios, not descent rules.
      - FM-SEC166 (main agent; `fm39/sec166_flip_descent_census.cpp`,
        `fm39/sec166_flip_single.cpp`).  FLIP DESCENT, a second exact
        descent next to (D).
        - For factors `u, v` of `Lambda` with `C = Lambda - u - v` and `A_ab(C)`
          the character table of `C`: `Q_uv = (Phi(Lambda) +
          Phi(Lambda^(uv)))/2 = sum_(c in CG(n_u,n_v)) Phi(C, eps_u eps_v c)`
          (`c = 0` counts `(1 + eps_u eps_v) Phi(C)`), and `Phi(Lambda) -
          Phi(Lambda^(uv)) = 4 eps_v A_(n_u n_v)(C)`.  So if `eps_v
          A_(n_u n_v)(C) >= 0` for some pair, `Phi(Lambda) >= Q_uv >= 0` by
          induction on the number of factors.  This is FM-MECH143's block
          identity at block size 2, used as a descent; it fails exactly at
          local minima of `Phi` under double sign flips.
        - Census (exhaustive): `Lambda = B + (sigma p)`, `sigma = (-1)^(minus
          count of B)`, so `FM3(Lambda) = 2 g_p(B)`.  Of the 21,364,508
          pair-free residual `(B,p)` through `W = 40` (994,058 cases where
          `B` holds `-sigma p` are excluded; they contain a pair), only 5,430
          have no flip descent.  By weight: 2, 4, 30, 14, 34, 33, 82, 66,
          146, 128, 274, 267, 416, 408, 686, 629, 1102, 1109 at
          `W = 23..40`.
        - On all 5,430, (D) holds, and the fixed rule TOPPAIR is monotone:
          remove the two factors of `B` of equal parity with the largest
          total label (ties: larger max label).  TopEven (the largest even
          label) is monotone whenever `B` has an even label; 8 no-flip lists
          have none, e.g. `B = (-1,3,5,7,9,11)`, `sigma p = -14`.
        - All eight large selector killers have a flip descent: W = 34
          (Rule M), 128, 131 (Rule 1'), 427, 540 and its reflection, 869
          (Rule W), and 1283 (Lemma WM), e.g. `(-11,+3)` at W = 1283.
        - CORRECTION 08:55 (found by luna_max_eris in FM-CHK96): for W = 131,
          427 and 1283, `B` has an odd number of minus signs, so the
          distinguished factor is `-p`.  The main agent had run
          `fm39/sec166_flip_single.cpp` with `+p` there.  Rerun with the
          correct sign, all three still have a flip descent, through
          different pairs:
          - W = 131, `sigma p = -15`: `(-1,-2)` = 18,354,066,022; the earlier
            `(-1,-7)` is negative.
          - W = 427, `sigma p = -31`: `(-1,-11)`.
          - W = 1283, `sigma p = -55`: `(-11,+3)`, with a different value
            from the earlier wrong-sign run.
          W = 34, 128, 540 (both) and 869 have an even minus count and were
          run correctly.  Other flip checks took `sigma p` from census logs
          or by parity and are unaffected.
        - CONJECTURE (FT): every pair-free residual `(B,p)` has a flip
          descent or a monotone TopPair removal.  With FM-MECH102,
          FM-MECH143 and the proved strata, (FT) gives the full FM3 cone.
          Neither part works alone: flip descent fails at sign minima
          (FM-MECH143's witness), and TopPair alone misses 261,978 residual
          cases at `W <= 40` (FM-SEC163).  Proof attempt: FM-MECH150
          (astra_max_vulcan); falsification: FM-SEC167 (luna_max_venus).
      - FM-SEC166 addendum (main agent; `fm39/sec166_regime_split_census.cpp`,
        all residual `(B,p)` through `W = 40`).  The two halves of (FT) live
        in separate regimes.  Let `w_TP` be the TopPair weight and `delta`
        the distance.
        - TopPair fails on only 126 cases, all with `w_TP <= 0.429 delta`,
          max label `<= 0.286 delta` and at least 15 factors.  On the 5,430
          no-flip lists, `w_TP >= 0.714 delta`, max label `>= 0.385 delta`,
          and there are at most 12 factors.
        - On the no-flip lists TopPair has a large margin:
          `g_p(B - TopPair) <= 0.0564 g_p(B)` (median 0.0159).
        - Split form of (FT), for a constant `c` in the gap (for example
          `c = 1/2`):
          - (T*) `w_TP >= c delta` implies TopPair is monotone;
          - (F*) `w_TP < c delta` implies a flip descent.
          Both hold on the whole census.  The gap narrows slowly with
          weight: the no-flip minimum of `w_TP/delta` falls from 1.25
          (`W = 23`) to 0.71 to 0.80 (`W = 38..40`), while the TopPair-failure
          maximum stays at 0.375 to 0.429.  Larger-weight data is needed.
        - CORRECTION 06:50 (found by astra_max_juno in FM-MECH151, confirmed
          by the main agent with the GMP evaluator `fm39/sec162_gp_single.cpp`):
          the two "kills" in the next bullet are WITHDRAWN.  The random search
          used signed 128-bit integers without an overflow guard, and those
          parents need 129..132 bits.  Exactly: TopPair is monotone on all
          four printed lists (`W = 124, 125, 141, 144`), with child/parent
          ratios 0.0000 to 0.0043.  The ratio split (T*)/(F*) and the
          few-ones split are therefore NOT refuted.  The census numbers are
          unaffected: their entries are at most `4^44 * 2(p+1)` (DP) and
          `2^N (T/N + 1)^(2N)` (Walsh), far below `2^127`.
        - Exact redo (GMP, `fm39/sec166_toppair_random_search_gmp.cpp`,
          W = 40..200, about 520,000 random residual lists):
          - The ratio split survives: every TopPair failure has `w_TP <=
            0.267 delta` (census: `<= 0.429 delta`), against `>= 0.714 delta`
            for every no-flip list at `W <= 40`.
          - The few-ones split (T1) is FALSE, exactly: `B = (-1)(2)^14
            (-3)^14 (-4)^2`, `sigma p = -11` (`W = 79`, 32 factors).  TopPair
            `(-4,-4)` raises `g_11` from 27,248,166,782,859,060,970 to
            34,145,565,434,211,036,460.  59 such failures have exactly one
            label 1.
          - Every failure checked has a flip descent (e.g. `(-1,-4)` here),
            so (FT) survives.
        - (Withdrawn; overflow) KILLED at 05:45 (main agent random search,
          `fm39/sec166_toppair_random_search.cpp` and `..._few_ones.cpp`):
          - (T*) with `c = 1/2` is false.  TopPair increases `g_p` at `W = 124`,
            `B = (+25,+7, 1's and 2's; 66 factors)`, `sigma p = 54`, with
            `w_TP/delta = 0.914`.
          - The split "at most one label 1 implies TopPair" is false at
            `W = 144` (`B = (-14, one 1, 2's and 3's; 56 factors)`,
            `sigma p = -44`).
          - Every such failure has a flip descent (e.g. `(+25,-1)`,
            `(-14,+2)`), so (FT) survives.  All TopPair failures found have
            many factors: at least 15 at `W <= 40`, and at least 26 in the
            random search over `W = 40..160`.  The open question is a split
            by factor count.
        - (FT) exhaustive through `W = 44` (main agent,
          `fm39/sec166_flip_descent_census_w48.cpp`, int128 safe since table
          entries are at most `4^44 * 2(p+1)`).  Of 66,116,748 residual cases,
          14,108 have no flip descent; on all of them (D) holds and TopPair
          is monotone.  No failure of (FT).
        - (FT) exhaustive through `W = 48` (main agent, same code, 25 min).
          Of 190,155,061 residual cases, 33,487 have no flip descent; on all of
          them (D) holds and TopPair is monotone.  No failure of (FT).
          - Structure of the no-flip lists:
            - factor counts 6..13 and 15 (6 lists with 15 factors, at
              `W = 45` and 48);
            - label-1 count 0, 1 or 2 (6 lists with two 1's);
            - label multiplicity at most 4.
          - KILL of lemma (F1) "two 1's give a flip" (GMP-checked by
            `fm39/sec166_flip_single.cpp` and the class-pattern tool
            `fm39/sec166_class_pattern_flips.cpp`): `B = (+1,+1,-2,+3^4,-4^3,
            +5^4)`, `sigma p = +8` (`W = 48`, `delta = 20`, 15 factors) has
            all 19 class-pair flips strictly negative.  TopPair `(5,5)` works
            (`g = 4,075,371`, child 160,741).
          - The no-flip minimum of `w_TP/delta` keeps falling: 0.588
            (`W = 41`), 0.556 (`W = 45`), 0.500 (`W = 48`, at this witness).
            TopPair failures reach 0.429 at `W <= 40`.  So the ratio split at
            `c = 1/2` sits at its edge and will likely break at larger
            weight.  (FT) itself is unaffected.
        - (FT) exhaustive through `W = 52` (main agent,
          `fm39/sec166_flip_census_resumable.cpp`: chunked, each finished chunk
          fsynced and skipped on restart, kill-and-restart test passed; 99
          min).  Of 516,000,611 residual cases, 75,532 have no flip descent; on
          all of them (D) holds and TopPair is monotone.  No (FT) failure.
          - Structure of the no-flip lists: factor counts 6..16 (13,282 with
            `>= 10` factors); at most two labels 1; multiplicity at most 4.
          - Minimum `w_TP/delta` = 5/11 (`W = 50`, the (F1) family).
        - Few-factor Walsh census (main agent,
          `fm39/sec166_walsh_few_factor_census.cpp`).  For each label
          multiset, one Walsh transform of `m(S) m(S^c)` gives `Phi` at every
          sign pattern, so all flip tests are table lookups.  It reproduces
          the DP census exactly at `W <= 40` (no-flip counts 499, 1597, 1944,
          1066, 270 for 6..10 factors).  Total weight `T = W + p`; no (FT)
          failure anywhere:

          | factors | `T <=` | no-flip patterns | TopPair failures |
          |---:|---:|---:|---:|
          | 6 | 250 | 7,707,257 | 0 |
          | 7 | 160 | 8,565,786 | 0 |
          | 8 | 124 | 5,471,268 | 0 |
          | 9 | 104 | 1,463,154 | 0 |
          | 10 | 92 | 247,753 | 0 |
          | 11 | 84 | 36,284 | 0 |
          | 12 | 80 | 6,380 | 0 |

        - TopPair on ALL pair-free patterns, not only no-flip ones
          (`fm39/sec166_walsh_toppair_all_patterns.cpp`): monotone for every
          list with 7..11 factors at `W <= 40`.  So at these weights TopPair
          fails only from 15 factors on, while no-flip lists have at most 12.
        - TopPair is monotone on EVERY residual pattern with 7 factors up
          to `T = 120` (81,513,734 patterns) and with 8 factors up to
          `T = 100` (121,763,604 patterns), not only on no-flip ones.
        - Structure of no-flip patterns (`fm39/sec166_walsh_noflip_structure.cpp`;
          757,869 with 8 factors at `T <= 100`, 46,037 with 10 factors at
          `T <= 80`): never two or more labels 1 (0 or 1 always); label
          multiplicity at most 4 (mostly 1 or 2); odd-label count 0, 4, 6, 8
          or 10, never 2.  Candidate lemma (F1): a residual list with two
          factors of label 1 always has a flip descent.
      - FM-SEC167 (luna_max_venus; `fm39/sec167_FT_local_search.py`).  No (FT)
        counterexample in 23 local-search starts at `W = 41..300` (every
        terminal still had a flip descent).  Weak evidence: the main-agent
        rerun did not finish, because one evaluator call exceeded the
        script's own 60 s subprocess timeout on the loaded machine.
        - Single fixed flips fail often.  For example, the flip of two equal
          labels `r` fails on 10% (`r = 1`) to 39% (`r = 4`) of lists with a
          repeated `r` at `W <= 36`.  So (F*) is an existence statement.
      - FM-MECH147 (astra_max_juno; `fm39/mech147_repeated_class_thresholds_repro.py`,
        rerun: ALL CHECKS PASS).  Exact terminating thresholds for
        repeated-class removal.
        - Exact numerator theorem for two `+1` or two `-2`: a translated
          polynomial test whose success proves FM3 and the pair removal for
          every larger multiplicity.  It cuts the FM-MECH144 thresholds
          from 460..1112 down to 2 or 3 on the listed cores (8 for
          `(-3)^64`).  Two `-2`: removal and FM3 hold for every multiplicity
          `>= 2` on three families.
        - Moment-ratio lemma: one successful exact test proves every larger
          count.  Implemented for `+-3` (`C = (-4,-4)`, `p = 6`: every even
          multiplicity `>= 494`) and `+2` (`C = (-3,+5)`, `p = 6`: `>= 266`).
          It gives effective eventual thresholds for every fixed signed
          label, with no uniform size bound.
        - Kill (knob: a cutoff independent of the remainder): no constant
          fundamental cutoff works (constant 4 is false on a residual word;
          exact endpoint limit `4(a^2-1)/25`).  A cutoff linear in the
          remainder weight is open: the first unresolved formula asks for
          `a >= c(p + |C| + 1)`.
      - FM-MECH148 (astra_max_ceres; `fm39/mech148_four_odd_dominance_repro.py`,
        rerun: PASS, identical output).  Four-odd layer.
        - For four odd labels on an even background `C`: `Phi = 2[H +
          e3 e4 J_1 + e2 e4 J_2 + e2 e3 J_3]`, with `H` the all-`x` term and
          `J_i` the three pairings.  The four values modulo odd reflection
          are `2(H +- J_1 +- J_2 +- J_3)` with an even number of minus signs.
          For the quartet (1,3,5,7) each of the three two-minus
          arrangements is the unique global minimum for some even
          background, e.g. (4,6,6), (2,6,8), (4,4,6).  So no minimizing
          prescription from the order of the odd labels alone exists (knob:
          background-independent sign order); 15 of 97 residual multisets
          have no minimum with all even signs positive.
        - Proved (Theorem 2): with `Psi_ij = U_i(x)U_j(y) - U_j(x)U_i(y)`,
          `Phi(D_P D_Q prod_A S prod_B S S_2^t) >= 0` whenever `P >= sum A`,
          `Q >= sum B` (same parity), for every `t`, nondecreasing in `t`.
          For odd `a, b` in `A`, the parent is at least `(min(a,b)+1)` times
          the child without `a, b`, with a nonnegative integer remainder.
          With two-odd fusion this is a positive descent; it covers
          residual sectors at levels 1 to 3 (e.g. `phi_3(h_7^2 h_2 h_4 h_6
          hat S_2^t) >= 2 phi_2(h_7^2 h_2 h_4 hat S_2^t)` for every `t`).
        - Kill (knob: removal of two odd labels only, including convex
          averages followed by fusion): at `B = (-1)^2 (+3)^2 (+2)^69`,
          `p = 18` (`W = 146`, `delta = 64`) every odd-pair removal increases
          `g_18`; removing one `+2` works.  (D) and FM3 survive.
        - Open: `H + min(J_1+J_2+J_3, J_1-J_2-J_3, -J_1+J_2-J_3, -J_1-J_2+J_3)
          >= 0` outside the dominance condition.
      - FM-MECH146 (luna_max_pluto; `fm39/mech146_channel_replacement_screen.py`,
        rerun: same output).  Replacement moves.
        - Same-sign channel replacement: replace two same-sign factors by
          one `+k`, `k in CG(n,m)`, `k >= 1`, if `g_p` does not increase.
          It is a valid descent and is implied by flip descent at that pair
          (all channel children are `>= 0` by induction).  Exhaustive
          through `W = 22`: 4,724 residual cases, no failure; least
          relative drop `34002/51874` at `B = (+1)^8 (-2)^4 (+3)^2`, `p = 6`.
        - Kills (knob: fixed selector): the smallest-duplicate selector at
          `B = (-11)^3 (-2)^2 (+6)^8`, `p = 11` (444,227,708 -> 549,660,463;
          `(-11,-11) -> S_2` works).  The same-sign correction-sign
          shortcut fails at `B = (-1,-2,-3^2,-4^2,-5)`, `p = 6`.  Main-agent
          note: that list still has a flip descent, but only through the
          pair `(-4, -6)` with the distinguished `-6` (difference 2).  So
          pairs that use `p` and opposite-sign pairs are needed.
      - FM-CHK95 (luna_max_eris, fresh code; `fm39/chk95_D_census_random_screen.py`,
        rerun: identical).  ACCEPT FM-SEC164 through `W = 36`: 1,530,800
        residual backgrounds, 6,590,648 `(B,p)` (this count includes the
        cases where `B` holds `-sigma p`), 59,919,062 removal comparisons,
        no (D) failure, all parents positive, minimum margin `22449/25937`.
        A random screen at `W = 100..400` (2,000 backgrounds, 99,533
        `(B,p)`, at most 12 factors, so only the few-factor regime) found no
        failure of (D), Rule W, Rule M or the most-frequent-pair rule.
      - FM-MECH149 (astra_max_minerva; `fm39/mech149_general_block_repro.py`,
        rerun with `--census`: PASS, table reproduced).  General even-block
        elimination.
        - Exact identity for any even block `I` (`|I| = 2k >= 4`, may contain
          `p`): all splits are bought by two-factor children with an integer
          correction `h_c`; `h_c >= 0` gives a positive descent.  Sufficient
          separation bound for `r = 2k >= 6`; knob: no bound `q >= f(D, M)`
          independent of `k` can make the correction coefficientwise
          nonnegative (`Q_11/Q_00 -> -4`).
        - Using the actual pure coefficients of `C` proves an infinite
          all-minus family (degree-10 polynomial identities).
        - Obstruction (knob: paying for every degree-permitted mixed channel,
          even adaptively): at one `W = 22` residual word all 55 proper even
          sub-blocks fail both tests; a mixed-sign rescue keeping the
          favorable mixed coefficient resolves it.
        - Coverage through `W = 30` (fully pair-free denominator 718,951):
          two-odd fusion, quartet separation and block tests reduce 89.37%;
          the uniform separation theorem adds only 0.137%.  For comparison,
          flip descent (FM-SEC166) reduces all but 265 of these cases
          (99.96%).
      - FM-MECH150 (astra_max_vulcan; `fm39/mech150_six_odd_FT_repro.py`,
        rerun: PASS).  A uniform (FT) sector, proved: background of six odd
        labels `a < b < c < d <= e <= f` (four smallest distinct), `p >= f`
        even, `f <= delta`.  Then TopPair removes `e, f` and `g_p(B) -
        g_p(B - TopPair) >= 4`, for every signing, all labels and
        distances.  It covers all eight census no-flip lists without an even
        label.
        - Exact graph reduction for six odd background labels: `g_p(B) = M +
          Q_A + s Q_W` (subsets of sizes 0, 2, 4 only), with explicit
          formulas for all 21 flip quantities.  Summing the 15 background
          flips gives `8(g_p - M)`, so a no-flip list there has `g_p < M`.
        - Normalization: `sum_(u<v) (Phi - Phi^(uv)) = 2 sum_S chi(S) |S||S^c|
          m(S) m(S^c)`.  This gives flip descent for pair-free lists with 2..5
          factors.
        - Knob for extending: unconditional TopPair monotonicity outside the
          pattern, e.g. `B = (-1)^4 (-2)(-3)^8 (-4)`, `p = 6`: TopPair
          `(-2,-4)` increases `g_6` from 996,550 to 1,077,706 (flip
          `(-1,-1)` works there).  Open: repeated labels among the four
          smallest, using the no-flip hypotheses.
      - FM-MECH153 (astra_max_minerva; `fm39/mech153_seven_eight_factors_repro.py`,
        rerun with `--receipts --census-log`: PASS, 142 receipt hashes
        match).  SEVEN FACTORS CLOSE; eight-factor 4|4 layer controlled.
        - Theorem 2 (new): every pair-free seven-factor word with all labels
          `>= 2` is nonnegative.  The largest negative 3|4 split `d` is paid
          by `N_7 + P >= 20 d >= T_3^-`, using interval bounds for the
          quartet's low channels (knob: ordinary fusion, no upper fusion
          boundary).
        - Theorem 3: every seven-factor word is nonnegative.  Words with a
          label 1 come from Corollary 23A9ZZ10 (two minus positions, every
          level) and Corollary 5A7B51 (the corrected shallow seven-factor
          residual) of `ginibre_q3/CENTRAL_CHARACTER_Q3_SEARCH.md`; the
          collector replays `PASS_EXACT_UNION` (308 tasks).
        - Eight factors: `N_8 >= T_4^-` (the whole negative 4|4 layer) for
          every pair-free word, so every eight-factor word with no negative
          3|5 term is nonnegative, in particular every all-odd one.  Knob:
          discarding the reserve `N_8 - T_4^-` fails at the no-flip list
          `(-1,+2,+3,+4,-5,+6,+7,+8)` by 9 (reserve 487).
        - Census: seven factors settle 1,597 of the 5,430 no-flip lists at
          `W <= 40`, and with six factors 2,096 (38.6%).
      - FM-MECH151 (astra_max_juno; `fm39/mech151_toppair_regions_repro.py`,
        rerun with the census audit: ALL CHECKS PASS).  TopPair regions.
        - Found that the W = 124 and 144 "TopPair failures" were signed-128
          overflow (see the CORRECTION in FM-SEC166).
        - Proved sufficient conditions for any same-parity removal:
          - the correlated class-parity condition (3), from unsigned fusion
            rows of submultisets;
          - an explicit few-factor threshold (4) in the minimum-label slack
            `h`, with an infinite nine-factor family;
          - every removal when `N <= 4`;
          - the whole boundary `w_TP = 2 delta`.
        - Census audit: the union covers 5,280 of 5,430 no-flip lists and
          none of the 126 TopPair failures.  The 150 left need signed block
          structure.  A stronger absolute-block inequality passed all 5,430.
      - FM-MECH152 (astra_max_ceres; `fm39/mech152_selected_flips_repro.py`,
        rerun: PASS).  Selected flips, proved:
        - Theorem 1, many fundamentals: `Lambda = (+1)^a C` with `a >= max(4r
          + 4, 2 kappa - 2r - 4)` (`kappa = (2r+3) lambda`, `lambda` quadratic
          in the core labels).  Flipping two `+1` is a descent, persisting
          for larger `a`.  This is a quantitative case of (F1).
        - Theorem 2: in `(-P)(-Q)(+q) prod(+h) (+2)^b` with `P >= sum h`, `q <
          Q`, the flip `(+q,-Q)` is a descent for every `b` (ordered cone,
          gap kept).
        - Kills (knob: prescribed label-weighted averages of fundamental
          flips, weights `n`, `n(n+2)`, `n-1`, `n^2-1`) at `W = 140` and
          `W = 393`.  Each witness still has a flip.
        - Census: no-flip lists have at most 12 factors and at most one
          label 1 through `W = 44`.  None lies in these regions, so the
          theorems do not shrink the census remainder.
      - FM-SEC168 (luna_max_pluto; `fm39/sec168_flip_pair_families.py`,
        `fm39/sec168_noflip_channel_replacement.py`; both rerun by the main
        agent with identical output, after fixing doubled regex backslashes
        in the extracted copy).
        - Every simple flip-pair family misses some flip-positive list in a
          22-list exact screen (W = 40..200): pairs with `sigma p`, with a
          smallest label, two copies of one label, opposite signs, (largest,
          smallest), maximum gap.  Minimal covering unions have two
          families, e.g. {pairs with `sigma p`, opposite-sign pairs}.  Small
          screen; not a census.
        - On ALL 1,904 no-flip lists at `W <= 36`, a same-sign channel
          replacement (FM-MECH146) strictly lowers `g_p`; minimum relative
          drop `352/429` at `B = (1,-2,3,3,-4,9,-10)`, `p = 12`
          (`(-4,-10) -> S_6`).  So channel replacement also covers the
          no-flip core there, next to TopPair.
      - FM-SEC165 (luna_max_mars; `fm39/sec165_D_census_w44_argmax.py`).
        - (D) exhaustive through `W = 44`: 68,797,116 residual `(B,p)` (this
          count includes the cases with `-sigma p` in `B`), 742,040,408
          removal comparisons, no failure, all parents positive.
        - Main-agent rerun capped at `W = 40` (the `W = 44` run needs about
          two hours and is not resumable): 22,358,566 cases and 222,096,442
          comparisons, matching FM-SEC164 exactly; no (D) failure.
        - Argmax catalogue: the most frequent best removal is the maximum
          label with its next same-parity label (45.5M of 70.3M argmax
          terms).  Rule W and Rule 1' each come within half of the best drop
          in every case through `W = 44`.  On the eight large witnesses, each
          fixed selector fails somewhere, but one of Rule W, Rule M, Rule 1'
          always works.
      - FM-SEC169 — WITHDRAWN 07:40 (main agent).  Its evaluator `entry`,
        and the "independent" `moment_A`, move a `+n` factor only on `x` and
        a `-n` factor only on `y`.  That is one assignment, not the
        character table of `prod (U_n(x) + eps U_n(y))`.  Check: for `B =
        (2..20)`, `p = 169`, the flip `(2,3)` is 243,634,987 by
        `fm39/sec166_flip_single.cpp` logic and by an unpruned GMP full table;
        FM-SEC169 printed 0.  Every value below is invalid; the search
        neither supports nor refutes (FT).  (Original text kept for the
        record.)  FM-SEC169 (luna_max_venus; `fm39/sec169_many_factor_noflip_search.py`,
        rerun: same output).  No no-flip residual list with many factors was
        found: 112 exact sign assignments on seven repeated-label profiles
        (20..80 factors, `W = 56..360`, multiplicity up to 30) plus random
        local searches all have a flip descent.  Many of the working flips
        have value exactly 0 (pair-creating flips on repeated labels;
        preserving flips on distinct-label runs).  Example with a strictly
        positive flip: `B = (+2)^8 (-3)^8 (-5)^8 (-8)^10 (-12)^11`,
        `sigma p = -268`, flip `(+2,+2)` = 26,942,364.  Finite sample; not
        a census.
      - TRIPLE-FLIP IDENTITY (main agent; `fm39/sec166_triple_flip_identity_check.py`,
        198 random cases, any `p`, signs and same-parity pair).  For a
        same-parity pair `(a, b)` of `B`, `C = B - a - b`, `Lambda = B +
        (sigma p)`, and the flip values `D_uv = eps_v A_(n_u n_v)(Lambda - u
        - v)`:
          `g_p(B) = T_1 + (D_ab + D_ap + D_bp)/2`,
          `T_1 = sum_(c in CG(a,b)) sum_(s in CG(p,c)) g_s(C)`.
        - Proof: of the four channel blocks of `g_p(B)`, the both-in-`x`
          block is `T_1`.  The other three blocks `T_2, T_3, T_4` satisfy
          `D_ab = T_3 + T_4`, `D_ap = T_2 + T_4`, `D_bp = T_2 + T_3`, using
          `G_C(t,s) = sigma_C G_C(s,t)`.
        - For a single even label `e`: `g_p(B) = sum_(s in CG(p,e))
          g_s(C) + D_(e,p)`.
        - Since `s = p` occurs in every `CG(p,c)` (`c` even, `c <= 2p`):
          `g_p(B) - g_p(B - R) = [T_1 - g_p(C)] + (D_ab + D_ap + D_bp)/2`,
          with `T_1 - g_p(C) >= min(a,b) g_p(C) + (other g_s(C)) >= 0` by
          induction.
        - Consequences:
          - a nonnegative triple sum `D_ab + D_ap + D_bp` makes that removal
            monotone, and `D_(e,p) >= 0` makes the removal of `e` monotone;
          - on a no-flip list, TopPair holds iff the same-side surplus
            `T_1 - g_p(C)` pays `(|D_ab| + |D_ap| + |D_bp|)/2`.
          So (FT) needs upper bounds on these three flips, not their signs.
      - FM-MECH157 (astra_max_ceres; `fm39/mech157_ratio_split_kill_repro.py`,
        rerun: PASS; witness cross-checked by the main agent with
        `fm39/sec166_flip_single.cpp` and `sec162_gp_single.cpp`).  KILL of
        the ratio split (knob: a fixed cutoff in `w_TP/delta` alone).
        - (F*) at `c = 1/2` is false: `Lambda = (-40,+42,-44,+46,+48,+50,
          +52,+54,+56,+58,+60,+62)`, `p = 62`, `W = 550`, `delta = 244`, `w_TP =
          118 < 122`.  All 66 flips are strictly negative (max
          -34,349,665).  TopPair works: `g_62 = 453,207,534,222,864`, child
          179,646,349,605.
        - Infinite no-flip family `n_i = 2(t+i)`, minus at indices 1, 2, `t >=
          1024`: `w_TP/delta -> 2/5`.  All flips are negative and TopPair is
          monotone there (degree-8 and degree-10 Newton certificates).
        - With the census TopPair failure at ratio `3/7` (`B = (-1)^13
          (-2)(-3)^5(-4)`, `p = 6`), no constant `c` splits the residual.
          (FT) survives.
      - FM-CHK97 (luna_max_mars, fresh code; `fm39/chk97_walsh_census_check.cpp`,
        main-agent rerun of `N = 6, T <= 120` and `N = 8, T <= 100/84`:
        identical).  ACCEPT the few-factor Walsh census: every requested count
        matches (residual multisets, patterns, no-flip patterns, zero TopPair
        failures, the all-pattern TopPair counts, the `W <= 40` no-flip
        counts).  Multiplicities come from inclusion-exclusion, and 80 sampled
        no-flip patterns pass a direct two-variable evaluation.
      - FM-MECH154 (astra_max_vulcan; `fm39/mech154_F1_sectors_repro.py`,
        rerun: PASS).  (F1) is false (the `W = 48` census witness).  What
        stands:
        - (F1) for at most six factors, and for a seven-factor sector (six
          odd background labels with two 1's, even maximum `p`);
        - (F_m) for at most six factors;
        - the shifted-label recurrence `(T+4) D_11 = 4(r-2) Phi((+1)^(r-2) C)
          + 4 sum_v (n_v+1) eps_v A_(1,n_v-1)((+1)^(r-2)(C-v))`;
        - a witness where every flip incident to a 1 is negative but
          `(+2,+2)` works, and one with multiplicity 5 where every
          equal-label flip is negative.
      - TopPair margin on no-flip lists through `W = 48` (main agent,
        `fm39/sec166_toppair_margin_noflip.py`, GMP evaluator, all 33,487
        lists): the worst `g_p(B - TopPair)/g_p(B)` per weight is at most
        0.0564 (`W = 27`) and at most 0.0414 for every `W >= 34`.  It does not
        grow with weight: on no-flip lists TopPair cuts `g_p` by at least
        94%.
      - FM-MECH155 (astra_max_minerva; `fm39/mech155_eight_factor_sectors_repro.py`,
        rerun with the census log: PASS).  Eight factors, three uniform
        sectors:
        - four odd labels, all labels `>= 2`;
        - six odd labels, all labels `>= 2`;
        - all-even labels, all labels `>= 4`.
        - Together with FM-MECH153 and two-odd fusion, every eight-factor
          word is covered if all labels are `>= 3`, or if all are `>= 2` and
          one is odd.  The proof uses five-factor channel certificates
          (722,142 endpoint cases, with an exact reduction of unbounded
          intervals) plus four finite exceptions.
        - Census: 306 of the 1,944 eight-factor no-flip lists, so 2,402 of
          5,430 in all.  The rest contain a label 1 (or are all-even with a
          2).  Open: `T_3^- <= (N_8 - T_4^-) + P + T_3^+ + T_4^+` there.
      - FM-SEC170 (luna_max_pluto; `fm39/sec170_ratio_split_screen.py`,
        rerun: same output).  Exact screens found no ratio-split failure
        (117,147 + 10,346 + 2,211 random residuals; 1,808 targeted high-ratio
        cases; 541 low-label flip checks).  Superseded by FM-MECH157, whose
        12-factor even-label witness refutes (F*).  Exact TopPair failures
        with flips: `W = 43` (ratio 1/3, `D_11 > 0`) and `W = 111` (ratio
        10/51).  Pluto also bounds the census arithmetic: table norms at most
        `2^125` through `W = 48`, so int128 is safe.
      - FM-CHK94 (luna_max_neptune; `fm39/chk94_D_falsification.cpp`; main
        agent reran `--mode check` and `--mode seeds` on the first seven
        witnesses with the same per-witness counts; the `W = 1283` full scan,
        about 4 h single-threaded, was not rerun).  No failure of (D) found:
        - Every one of the eight ledger witnesses has a monotone removal
          (1,037 removals scanned in full).  The removals are mostly
          monotone: 35/36 at `W = 128`, 129/130 at `W = 427`, 161/162 at
          `W = 540`, 135/136 at `W = 869`.  The smallest best margin is at
          `W = 34` (`476553/498275`).
        - Thirteen structured rows at `W = 500..3000` and eight one-step
          mutants of the witnesses also pass.  Finite evidence only.
      - FM-CHK98 (luna_max_mars, fresh code; `fm39/chk98_triple_identity_mech157_check.py`,
        `fm39/chk98_gmp_crosscheck.py`; both rerun by the main agent:
        PASS).  ACCEPT:
        - the triple-flip and single-even-label identities, with Mars's own
          proof and 3,000 random exact checks each (both minus parities,
          repeated labels, equal-label pairs, pairs with `p`'s label);
        - FM-MECH157's 12-factor witness and its `t >= 1024` family
          certificates;
        - the ratio-3/7 TopPair failure, which has a flip descent in every
          signing.
      - FM-SEC171 (luna_max_venus; `fm39/sec171_class_profile_FT_screen.py`,
        `fm39/sec171_F1_direct_check.py`; the direct check rerun by the main
        agent: same values).  With the corrected evaluator:
        - (F1) is refuted independently at the `W = 48` witness (all 19
          flips listed).
        - Exact class-pattern scan of 97 profiles (292 profile-`p` runs,
          17,944 patterns, `W = 48..200`) finds six no-flip patterns at
          `W = 48, 50, 56` (15 or 16 factors, `w_TP/delta` down to 5/11).
          TopPair is monotone on every one; none was found at `W >= 100`.
        - The FM-SEC169 zeros came from the withdrawn evaluator, e.g. the
          true `A_22((+2)^28 (+3)^30 (+4)^20 (+214)) = 118,588,718,722`.
      - FM-MECH156 (astra_max_juno; `fm39/mech156_moment_certificate_repro.py`,
        rerun: ALL CHECKS PASS).  A proved rational removal certificate (E).
        - Lemma 1: exact moments of the positive rotation kernel `M_c`
          (positive semidefinite for `0 <= c <= 1`), checked on 19,855
          identities.  Lemma 2: the four channel blocks of a same-parity
          removal come from one endpoint matrix.
        - Proposition 3: an explicit rational test (E) (Cauchy-Schwarz in the
          kernel's Hilbert space) implies `g_p(B) >= g_p(B - R)` for any
          same-parity pair `R`, uniformly in labels, signs, factor count and
          distance.  If `w_TP >= delta`, it simplifies (no wraparound term).
        - At `W <= 40`, (E) holds on all 5,430 no-flip lists (including the
          150 missed by FM-MECH151) and rejects all 126 + 64 TopPair failures.
        - Knob kills: discarding the favorable mixed block fails; so does
          monotonicity under growth of the removed labels.
        - Open: (T*) with `c = 1`, i.e. (E') `w_TP >= delta` implies TopPair.
        - Main-agent audit at `W <= 48` (`fm39/mech156_certificate_audit_w48.py`):
          (E) holds on 33,485 of the 33,487 no-flip lists (count corrected
          by Juno in FM-MECH161: four of the six failures are the W = 50/56
          lists the main agent added).  It FAILS on the six
          members of the (F1) family (`+-(1^2, 2, 3^4, 4^3, 5^4)`, `sigma p = 8`
          at `W = 48`; with one more `-2`, `p = 6` at `W = 50`; with `+8`, `p = 8`
          at `W = 56`), where TopPair is still monotone.  The 12-factor
          FM-MECH157 witness (`W = 550`) exceeds the script's size limit and
          was not tested.
      - FM-SEC173 (luna_max_pluto; `fm39/sec173_w550_noflip_check.py`, rerun:
        same values).  Exact screens at `W = 50..150` (61 multisets, every
        sign pattern; 192 local-descent starts) found no no-flip list and no
        TopPair failure.  An independent Python check confirms the W = 550
        FM-MECH157 witness (all 66 flips negative, TopPair monotone).
      - FM-MECH158 (astra_max_ceres; `fm39/mech158_translated_families_repro.py`,
        rerun with the W <= 48 census: PASS).
        - Proposition 1, a sufficient triple-flip bound: the same-side
          surplus `T - g_p(C)` dominates an unsigned bound on the three
          mixed blocks.  It implies TopPair monotonicity; no no-flip
          hypothesis is needed.
        - Theorem 2: TopPair monotone and FM3 on 350 uniform translated
          families with arbitrary pair-free signs:
          - `L` = 6..12 labels in arithmetic progression, step `q` in {1,2};
          - at most one label perturbed by `+-2`;
          - all translates `t >= 20` (`q = 1`) or `t >= 0` (`q = 2`).
          This includes the whole FM-MECH157 family.  Proof: polynomial
          certificates in `t`.
        - Census audit (33,487 no-flip lists, `W <= 48`):
          - Proposition 1 certifies 32,955;
          - the signed bound `T - g_p(C) >= |X| + |Y| + |Z|` holds on ALL
            33,487;
          - 532 are left by the unsigned test.  Knob: replacing the signed
            background by its unsigned table, first lost at `W = 29`.
      - FM-MECH159 (astra_max_vulcan; `fm39/mech159_mixed_certificate_repro.py`,
        rerun: PASS).  Mixed removal-plus-flip certificate.
        - With `C_R = Delta_R + (1/(2L(L-1))) sum_(u<v) D_uv`, both `C_R` and
          `Delta_R` are at least a positive multiple of `m(Lambda)` in these
          cases: `L = 7` with minimum label `>= 12`, `L = 8` with `>= 22`,
          `L = 9` with `>= 32`.  So every same-parity removal is monotone
          there.  FM3 holds for `L = 8`, minimum label `>= 20`, and for `L = 9`,
          minimum label `>= 28`.
        - Low-channel lemma: `mu_(2r)(A) >= (r+1) mu_0(A)` for 3 or 4
          factors, `2r <= min A`.  The 5-factor case (needed for `L = 10, 11`)
          passed 9,268 checks and is open.
        - Kill (knob: a constant-coefficient blend): every fixed blend of
          TopPair with the uniform global flip sum fails somewhere.
      - FM-SEC172 (luna_max_neptune; `fm39/sec172_run_family_FT_screen.cpp`;
        main-agent rerun of `--check-only` and `--next-numeric-only`: identical,
        9,857 profiles, 8 no-flip, 0 TopPair failures).  Runs of 8..20
        consecutive labels (step 1 or 2) with 2/4/6 minus signs at all
        positions: 31 no-flip lists found, TopPair monotone on all; the
        worst child/parent ratio is 3343/190903 = 0.0175, at
        `(-1,+2,+3,+4,-5,+6,+7,+8,+9,+10)`, `p = 11`.  Starts `s >= 20` are
        only sampled (time limits).
      - FM-SEC171 addendum: the main-agent rerun of the full 97-profile
        class-pattern scan finished: PASS.
      - FM-MECH160 (astra_max_minerva; `fm39/mech160_eight_factors_close_repro.py`,
        rerun: FM-MECH160 PASS, identical output; FM-CHK102 ACCEPT).
        EIGHT FACTORS CLOSE: with FM-MECH153/155 and two-odd fusion, FM3
        holds for every list with at most eight factors.
        - Remaining sectors: four or six odd labels with a label 1, and
          all-even labels with a label 2.
        - If the four smallest labels are all `<= 4`, the residual condition
          and the FM-MECH142 separated-quartet descent (children with at
          most six factors) leave a finite set: every label `<= 48`.  It is
          checked exactly: 14,903 unsigned words, 517,426 signings, minimum
          `Phi = 270` at `(1,1,1,-2,3,-8,-8,-8)`.
        - Otherwise at most three labels are `<= 4`, so at most one triple
          is all-small, and uniform triple-profile and fivefold
          row certificates pay the 3|5 layer.  Tiles: 1,189,029 + 1,165,441
          + 3,936,471 + 575,995 + 591,843 + 1,480,633, all passing.
        - Open, nine factors: `T_(3|6)^- + T_(4|5)^- <= N_9 + P_9 +
          T_(3|6)^+ + T_(4|5)^+`.
      - FM-MECH161 (astra_max_juno; `fm39/mech161_balanced_certificate_repro.py`,
        rerun: ALL CHECKS PASS).  Balanced certificate (CB), proved.
        - Lemma 1 keeps the wraparound term exactly: `X + Y = A - B_0`.
        - Proposition 2: for any `0 <= k <= 2b`, a tilted Cauchy-Schwarz bound
          `A^2 <= E_k` gives the rational test `Q >= 0`, `Q^2 >= min_k E_k`.
          It implies `g_p(B) >= g_p(B - R)` for any same-parity pair `R`,
          uniformly, with no no-flip hypothesis.  The `k = 0` case dominates
          (E).  Knobs that (E) loses: unbalanced norm allocation and
          charging the wraparound term with the adverse sign.
        - Audit: (CB) holds on ALL 33,487 no-flip lists at `W <= 48`, on the
          W = 50/56 (F1)-family lists and on the W = 550 FM-MECH157
          witness.  It rejects all 190 recorded TopPair failures.
        - Kill (knob: using only the three flips of the triple identity):
          at a `W = 34` list all three are negative and TopPair fails, while
          another flip is positive.
        - Open: no flip descent implies `Q >= 0` and `Q^2 >= min_k E_k`.
      - FM-CHK99 (luna_max_mars, fresh code; `fm39/chk99_seven_eight_factor_check.py`,
        main-agent rerun with `--seven-max 16 --word-max 12 --eight-max 9`:
        PASS).
        - FM-MECH153 Theorem 2: REPAIR of scope only.  The payment
          inequality `N_7 + P >= 20 d >= T_3^-` holds for even minus parity,
          which is all the proof uses (odd-parity words vanish).  For all
          signings it is false: `(-2)^7` has `N_7 = 36`, `P = 126`, `d = 3`,
          `T_3^- = 105`, `Phi = 0`.  Census: 2,340,495 even-minus signings
          with labels 2..16, no failure.
        - FM-MECH153 Theorem 3: ACCEPT.  The reductions land in Corollaries
          23A9ZZ10 and 5A7B51; the receipt replay passes.  All 2,035,800
          signed seven-label multisets with labels 1..12 have `Phi >= 0`
          (minimum 0 at `(1^6, 8)`).
        - FM-MECH155: ACCEPT (tiling, channel bounds, bounded census with
          labels `<= 9`).
        - FM-MECH160 (the full eight-factor closure) came later and is not
          covered by this check.
      - FM-CHK100 (luna_max_pluto, fresh code; `fm39/chk100_certificate_E_check.py`,
        main-agent rerun after fixing an indexing typo `mat[i,j]` ->
        `mat[i][j]` in the extracted copy: ALL EXACT CHECKS PASS).  ACCEPT
        FM-MECH156 Lemma 1 (kernel positivity, moment formula) and the
        implication (E) => monotone removal.  The induction fact used is
        only `g_p(C) >= 0`.  REPAIR of the count: (E) holds on 33,485 of
        33,487 at `W <= 48`.
      - FM-CHK96 (luna_max_eris, fresh code; `fm39/chk96_flip_census_check.py`;
        main-agent rerun: identical totals, per-W counts and ratio 31/550).
        - ACCEPT the FM-SEC166 census at `W <= 40`: 22,358,566 raw cases,
          994,058 exclusions, 21,364,508 residual; 5,430 without flip
          descent, with the stated per-W counts; (D) and TopPair hold on all
          of them; maximum TopPair ratio `31/550`.
        - REPAIR (scope): the formula `Phi - Phi^(uv) = 4 eps_v A` needs even
          total minus parity of `Lambda`, which every census list has.  In
          general use `A = (1 + sigma_C eps_u eps_v) T_ab / 2`.
        - REPAIR: the W = 131 pair `(-1,-7)` is not a descent (already
          corrected at 08:55); `(-2,-3)` also works.
      - FM-MECH162 (astra_max_ceres; `fm39/mech162_locality_obstructions_repro.py`,
        rerun: PASS).  Obstructions to proving (FT) from partial data.
        - Every flip is a signed cut sum of the subset expansion.  The three
          TopPair/`p` flips do not imply TopPair (exact residual example,
          drop -22,854,900).
        - Theorem 3: no criterion using flips among any fixed number of the
          largest factors implies TopPair.  Family `(+1)^A (-n_1)...(-n_2h)
          (+p)`, `A` large: all flips among the cores and `p` are negative and
          TopPair increases `g_p`, but `(+1,+1)` is a positive flip.  Exact
          for `A >= 1024` via 13 polynomial certificates.  Knob: locality of
          the flip hypotheses.
        - Proposition 4: flips never see `m(Lambda)`.  Lowering it from 624
          to 104 in an artificial table keeps all 28 flips negative, every
          proper even-sign value nonnegative and every full Walsh value
          positive, yet makes the TopPair drop negative (`S = -143`).  So
          flip signs plus positivity of shorter lists cannot imply (FT); the
          full-to-child Clebsch-Gordan relation must be used.  Knob:
          forgetting that relation.
        - Census: `S >= 0` on all 33,487 no-flip lists, with minimum
          `S/g = 377/550`.  At the extremes the nearest-to-zero flip lies
          outside TopPair and `p`.
      - FM-SEC174 (luna_max_venus; `fm39/sec174_run_family_margin_screen.py`).
        Runs of 8..16 labels (step 1 or 2, starts `s <= 200`, `W = 60..300`,
        2..6 minus positions, sampled `p`):
        - 967 no-flip sign patterns, all TopPair-monotone;
        - worst ratio `84983/2125083` (about 0.04), at the `W = 50` (F1)-family
          list, so no violation of (FT') "child <= parent/16";
        - maxima cross-checked with `sec166_flip_single.cpp` and
          `sec162_gp_single.cpp`.
        Finite search; `p` is sampled.  Main-agent rerun: all scan bands
        identical (88/319/420 and 52/74/14 no-flip, shape scan 3).  The
        script's final repository cross-check step crashed (TypeError,
        apparently a rendering defect in the printed code).
      - FM-MECH163 (astra_max_vulcan; `fm39/mech163_removal_thresholds_repro.py`,
        rerun with the census: PASS).
        - Theorem 1, low-channel bound: if at least three labels of `A` are
          `>= 2r`, then `mu_(2r)(A) >= (r+1) mu_0(A)`, with any number of
          factors.  This proves the five-factor lemma of FM-MECH159.  Knob:
          two large labels do not suffice (`A = (1,10,11)`, `r = 2`).
        - Unconditional removal descent (any same-parity removal monotone,
          given the shorter-list induction): for 7, 8, 9, 10, 11 factors
          when the minimum label is `>= 3, 7, 11, 15, 21`; in general the
          threshold is `O(2^(L/3))`.  These are descent sectors, not
          standalone strata.
        - Census: the label-based certificate covers 23,739 of the 33,487
          no-flip lists at `W <= 48`.  The rest have small minimum labels.
          Open: no flip descent implies the aggregate signed partition sum
          is `>= g_p(B - TopPair)`.
      - FM-CHK102 (luna_max_mars, fresh code; `fm39/chk102_eight_factor_closure_check.py`,
        rerun: same output).  ACCEPT FM-MECH160 on all three items:
        - the reduction to two sectors;
        - the cutoff `M <= D`, `q <= 2D`, `p <= 3D` (all labels `<= 48` when
          `D <= 16`), with separated-quartet children of at most six
          factors, and an independent finite census of 14,903 words and
          517,426 signings with minimum 270;
        - the uniform certificates, including the validity of the interval
          cap at 41.
        - Also exhaustive: all 6,435 label profiles `1..8` with all
          1,647,360 signings are nonnegative; 300 random words with labels
          up to 400 pass.
      - FM-CHK103 (luna_max_pluto, fresh code; `fm39/chk103_balanced_certificate_check.py`,
        rerun: ALL EXACT CHECKS PASS).  ACCEPT FM-MECH161:
        - Lemma 1 (`X + Y = A - B_0`);
        - Proposition 2 (tilted Cauchy-Schwarz; only `g_p(C) >= 0` is
          used);
        - the audit: all 33,487 no-flip lists certify, and all 190 TopPair
          failures fail (CB) for every `k`;
        - 638 random residual/pair cases, every (CB)-certified removal
          monotone.
      - FM-MECH166 (astra_max_ceres; `fm39/mech166_spread_label_region_repro.py`,
        rerun: PASS).  Uniform minimum-label region, no upper bound on the
        maximum label `M`.
        - Lemma 1, a coefficientwise low-channel bound: for three labels
          `>= m`, `U_(2r) U_a U_b U_c` dominates `c_m(r) U_a U_b U_c`.  Hence
          `mu_(2r)(A) >= (r+1) mu_0(A)` whenever `A` has three labels `>= 2r`,
          for any factor count.
        - Lemma 2 gives closed-form bounds for the proper-subset
          contributions.
        - Thresholds, with `L` counting the distinguished factor:

          | `L` | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
          |---|---|---|---|---|---|---|---|---|---|---|---|
          | FM3 directly | 3 | 4 | 6 | 8 | 10 | 13 | 17 | 22 | 28 | 35 | 45 |
          | every same-parity removal monotone | 4 | 4 | 6 | 8 | 11 | 14 | 18 | 22 | 28 | 36 | 45 |

        - Knob: bounding proper-subset contributions separately.  The
          refined certificate leaves 29,096 census no-flip lists (the
          covered ones have at most eight factors).  Open: the no-flip
          inequalities must recover the discarded cancellation.
      - FM-MECH164 (astra_max_minerva; `fm39/mech164_nine_factors_close_repro.py`,
        rerun: FM-MECH164 PASS, identical; FM-CHK106 ACCEPT).  NINE
        FACTORS CLOSE: FM3 holds for every list with at most nine factors.
        - Layer identity `Phi/2 = N_9 + P_9 + T_(3|6)^+ - T_(3|6)^- + T_(4|5)^+ -
          T_(4|5)^-`.
        - When at most four labels are `<= 6`, uniform five- and sixfold rows
          come from coefficient inequalities: `U_(2s) prod_A U_n >= H^(7)_s
          prod_A U_n` when three labels are `>= 7`, and `H^(8)` when all labels
          are even and three are `>= 8`.
        - With more small labels, the residual cutoff and separated-quartet
          descent reduce to a finite check: 3,049,184 patterns and
          160,310,765 signings, minimum `Phi = 578` at
          `(1,1,1,1,-2,-2,8,8,8)`.
        - Also: ten factors with minimum label `>= 16`.
        - Census: all 7,860 nine-factor no-flip lists are settled, so 28,911
          of 33,487 are proved outright.  4,576 remain, all with `>= 10`
          factors.
      - FM-MECH165 (astra_max_juno; `fm39/mech165_cb_band_repro.py`, rerun:
        ALL CHECKS PASS).
        - Lemma 1: every flip is an adjacent difference of kernel moments of
          its own deleted background.  So different flips use different
          kernels, and treating them as one matrix loses compatibility.
        - Theorem 5 (proved under shorter-word induction): every pair-free
          residual word with at least nine factors and `w_TP >= 2 delta - 2`
          satisfies (CB) with `k = 0`.  So TopPair is monotone there and
          FM3 holds, uniformly in labels, with no flip hypothesis.  This uses
          Lemma 3 (first two kernel diagonal rows) and Lemma 4 (Bernstein
          bound for the degree-4 diagonal).
        - Census: minimizing over every `k` in `0..2b`, all 33,487 pass; the
          minimum relative margin is about 0.45.  The band covers 1,024 of
          the 12,436 no-flip lists with at least nine factors.
        - Witnesses (knob: the set of flip hypotheses): neither the
          `p`-incident flips plus the TopPair flip, nor in addition all
          opposite-parity flips, imply (CB).  The full all-flips hypothesis
          is needed.
        - Open: `w_TP <= 2 delta - 4` with all flips negative implies (CB).
      - FM-CHK105 (luna_max_mars, fresh code; `fm39/chk105_thresholds_obstructions_check.py`,
        rerun: PASS).  ACCEPT, with no repair:
        - FM-MECH163 Theorem 1 and the removal-descent thresholds, with
          bounded-window descent checks;
        - FM-MECH162 Theorem 3 (locality obstruction) and Proposition 4
          (artificial-table obstruction).
      - FM-MECH167 (astra_max_vulcan; `fm39/mech167_small_flip_obstruction_repro.py`,
        rerun: ALL EXACT CHECKS PASS).  Obstruction, knob: support of the flip
        hypotheses.  Families `(+1)^r (+2)^A (-3)^4 (-4)^2 (+(8-r))` and
        `(+1)^r (-2)^2 (-3)^4 (+4)^A (-6)^2 (+(10-r))`, `r` in {1,2}.  In
        both:
        - every flip involving label 1 is negative, and in the second family
          every flip involving label 1 or 2;
        - the three TopPair/`p` flips are negative;
        - yet `3/2 < g_p(B - TP)/g_p(B) < 2`, so TopPair fails.
        A flip of two copies of the repeated class (`(+2,+2)`, `(+4,+4)`) is
        positive, so (FT) holds.  Exact witnesses at `A = 128, 192`; the
        signs hold for large `A` (corner asymptotics, ratio about
        `A^2/13500` and `A^2/27440`).  So no certificate supported on small
        flips plus the TopPair triangle can prove (FT).
      - FM-SEC175 (luna_max_venus; `fm39/sec175_all_p_run_mixed_screen.py`;
        main-agent rerun: identical, 18,418 + 22,954 no-flip patterns).  Every admissible `p`, `W = 60..200`:
        - runs of 8..12 labels, and runs plus `1^a, 2, 3`: 17,323 runs and
          20,156,928 sign patterns;
        - 41,372 no-flip patterns, all TopPair-monotone, every
          child/parent ratio `< 1/16`, so (FT') holds there;
        - maxima cross-checked with `sec166_flip_single.cpp` and
          `sec162_gp_single.cpp`.
      - FM-CHK101 (luna_max_neptune, fresh code; `fm39/chk101_mech158_mech159_check.cpp`,
        rerun: ALL CHECKS PASS).
        - ACCEPT FM-MECH158: the translated families, with their
          stabilization, degree bounds and Newton certificates.
        - ACCEPT FM-MECH159: the 7/8/9-factor thresholds.  This includes the
          exhaustive box at `L = 9`: all 5,005 multisets at the threshold,
          57,799 even signings, all positive.  The five-factor low-channel
          estimate is not part of the accepted claim; it was later proved
          as FM-MECH163 Theorem 1 and FM-MECH166 Lemma 1.
      - FM-MECH168 (astra_max_ceres; `fm39/mech168_label_ge2_region_repro.py`,
        rerun: PASS).
        - Lemma 1: `d(Lambda) >= Gamma_R d(Lambda')`, with `Gamma_R = sum_r
          H_(Lambda')(r) >= b + 1`.  It keeps all parent-to-child fusion
          channels via lower profiles `H_A(r)` (pairs, exact triple profile,
          third-largest-label profile from FM-MECH166 Lemma 1, quartet
          interval bound).
        - Theorem 2: for pair-free lists with every label `>= 2`, an
          explicit label-and-sign quantity `rho_R` satisfies `g_p(B) -
          g_p(B - R) >= rho_R d(Lambda)`.  So `rho_R >= 0` gives removal
          monotonicity and FM3 directly, with no shorter-list assumption.
          Corollary 3: an explicit region with at most two labels below
          `M` and arbitrarily many other factors.
        - Census (no label 1, minimum label 2 or 3, `>= 9` factors, `W <= 48`):
          798 of 832 covered, including all 736 with nine factors.  34
          ten-factor lists remain.
        - Kill (knob: locality): flips involving only the minimum labels
          (2 or 3, even repeated) do not imply TopPair; exact residual
          witnesses for both minima.
      - FM-MECH169 (astra_max_minerva; `fm39/mech169_ten_factors_close_repro.py`,
        rerun: FM-MECH169 PASS, identical counts, 221 s; FM-CHK108
        ACCEPT).  TEN FACTORS CLOSE: FM3 holds for every list with at most
        ten factors.
        - Payment identity with the 3|7, 4|6 and 5|5 layers.
        - Uniform low-channel profiles `H^(9)` and `H^(10)` (all-even), from a
          finite tiling catalog of 8,930 and 6,336 shapes, covering
          unbounded labels.
        - Lemma 3 keeps actual `U_2` factors for the all-even payment.
        - Finite box: 100,830,299 label patterns and 9,728,129,290 signings,
          all nonnegative.
        - Census: all 3,318 ten-factor no-flip lists at `W <= 48` are
          settled; 32,229 of 33,487 are proved outright, and 1,258 with `>= 11`
          factors remain.  Next: eleven factors (uniform payment with at most
          six labels `<= 8`).
      - FM-CHK104 (luna_max_eris, fresh code; `fm39/chk104_cb_audit_w52.py`,
        rerun: same output).  ACCEPT for the `W <= 52` census.
        - All 75,532 no-flip lists pass the balanced certificate (CB) for
          TopPair, so no other removal is needed.  TopPair is monotone on
          all of them.
        - 512 sampled lists were rechecked as no-flip with a fresh
          evaluator.
        - Worst integer (CB) margin 36, at `B = (-1,-3,-6,-7,-8)`, `sigma p =
          -9`.  Largest TopPair ratio 31/550.
        - So through `W = 52` every residual list has a flip descent or a
          (CB)-certified TopPair descent.  Open: the same for `W > 52`.
      - FM-CHK106 (luna_max_mars, fresh code; `fm39/chk106_nine_factor_closure_check.py`,
        rerun: PASS, the same finite census).  ACCEPT FM-MECH164 (nine
        factors).  Repair: the separated-quartet children have at most
        seven factors, covered by the proved results.  The independent
        finite census reproduces 3,049,184 patterns, 160,310,765 signings
        and the minimum 578.
      - FM-MECH171 (astra_max_vulcan; `fm39/mech171_multiplicity_threshold_repro.py`,
        rerun: ALL EXACT CHECKS PASS).
        - For `Lambda = (eps n)^r C` with `t` other factors, largest label
          `Q`, `q = Q + 1`: if `r >= 2 + 128 q^2 chi(n,eps)(t+3) ceil(log2 q)`
          (`chi = 1` unless `eps = -1` and `n` is even, then `(n+1)^2`), the flip
          of two copies of `eps n` is positive and `Phi(Lambda) > 0`.  This
          holds uniformly in the other factors.
        - Knob: the threshold depends on the remainder; an absolute
          cutoff (multiplicity `>= 5`) is open.  The successful flip need
          not involve a designated class of multiplicity five (exact
          obstruction).
      - REALLOCATION (user decision, 2026-10-02, 11:00).  Case-by-case
        grinding stops: factor-count closures after eleven factors, census
        extensions, and sector certificates from coefficient bounds.  The
        agents work on structure:
        - FM-STR1: an injection or positive-operator mechanism for the
          involution trace `Phi = dim Inv^+ - dim Inv^-`;
        - FM-STR2: a group or duplicate-variable change that makes the
          level-`r` integrand a genuine character;
        - FM-STR3 and FM-STR4: an exact identity "TopPair drop + sum
          alpha D = manifestly nonnegative", which would give the no-flip
          implication;
        - FM-STR5: exact computational support;
        - the bounded factor count of no-flip lists, and one checker.
      - FM-CHK107 (luna_max_neptune, fresh code; `fm39/chk107_mech166_mech165_check.py`,
        rerun: PASS).  ACCEPT, with no repair:
        - FM-MECH166: Lemma 1, the closed forms, the thresholds, boundary
          boxes at `L = 9, 10` and large-label samples;
        - FM-MECH165 Theorem 5: the band `w_TP >= 2 delta - 2`, `>= 9`
          factors, under the shorter-word assumption.
      - FM-STR3 (astra_max_juno; `fm39/str3_one_generation_obstruction.py`,
        rerun: ALL EXACT CHECKS PASS).  Structural obstruction, knob: depth of
        the fusion constraints.
        - Root: an actual eleven-factor residual list, `p = 13`, `W = 55`,
          `delta = 21`, TopPair `(+8,+10)`.  Lower `m(Lambda)` (222,024 ->
          30,265) and the invariants of its 240 one-pair fusion children by
          an exact integer certificate.
        - What stays the same: all 55 flips (max -185), every parent scalar
          and signed fusion identity (56,320 equations), nonnegativity of
          every shorter signed value, and all multiplicities.
        - Yet the TopPair difference becomes -758.
        - So flips + shorter-list positivity + ALL one-generation fusion
          relations do not imply TopPair.  The table violates a
          child-to-grandchild fusion relation (defect 68,353), so an
          identity must use recoupling consistency across common deleted
          backgrounds, at least two generations deep.  Whether two suffice is
          open.
        - Stopped FM-MECH170 (band widening): the band `w_TP >= 2 delta - 2`
          stands; the attempted extensions are unverified.
      - FM-STR2 (astra_max_vulcan; `fm39/str2_haar_model_obstruction.py`, rerun:
        PASS).  Structural no-go, knob: a common Haar pushforward with
        positive-definite realizations.
        - For every level `r >= 1` there are two actual pair-free word
          integrands `A, B` with `Cov_(nu_r)(A, B) < 0` under the tilted measure
          `nu_r = (x-y)^(2r) mu x mu / Z` (exact moments by integration by parts;
          e.g. `Cov = -1` at `r = 1`, `Cov(S_1^2, S_2) = -7/25` at `r = 2`).
          Real positive-definite functions on a compact group have
          nonnegative covariance under Haar.  So NO compact group of any type
          or rank, with any map pushing Haar to `nu_r`, makes all word
          integrands positive definite.  This includes genuine characters and
          nonnegative real combinations of characters, and spherical
          functions on homogeneous spaces.
        - A single literal character fails by integrality: for the pure-minus
          word `(-3)(-1)^(2r-1)`, its `nu_r`-mean is strictly between 0 and 2
          and never 1, for `r >= 2`.
        - Ginibre-type angles `x = 2cos(alpha+beta)`, `y = 2cos(alpha-beta)`
          make every `S_n` and `-D_n` a positive-semidefinite kernel.  But
          the universal contraction is not a squared norm feature by
          feature (`S_2` has feature contributions +2 and -2), and `(x-y)^(2r)`
          has no measure-only Gram factorization (determinant -1).  A Gram
          proof must use cancellation across the complete tensor product.
      - FM-MECH173 (astra_max_minerva; `fm39/mech173_eleven_factors_close_repro.py`,
        rerun with 24 threads: FM-MECH173 PASS, identical counts, 520 s;
        independent check pending).  ELEVEN FACTORS CLOSE: FM3 holds for
        every list with at most eleven factors.
        - The small-label cutoff moves from 8 to 12.
        - A new four-large-label coefficientwise profile `H_4^(13)` (`H_4^(14)`
          all-even) comes from a finite tiling reduction.
        - Most seven-small-label prefixes are certified uniformly.
        - Finite box: 326,484,586 label patterns and 47,094,362,366 signings,
          all nonnegative.
        - Census: 888 more no-flip lists at `W <= 48` are settled (370 with
          `>= 12` factors remain), and 2,544 at `W <= 52` (1,330 remain).
        - The factor-count route STOPS here by user decision (the finite boxes
          grow too fast to reach every factor count).
      - FM-STR1 (astra_max_ceres; `fm39/str1_block_injection_repro.py`, rerun:
        FM-STR1 PASS).  A structural mechanism on a narrow class, plus
        locality no-gos.
        - Theorem 1: for two minus factors `-n, -m` with plus weight `P <=
          n + m`, `I^- = 0` if `P < n + m`.  At `P = n + m` there is an explicit
          `G`-equivariant, `tau`-odd block-transfer operator `J` (singlet
          transfer of the block {`-m`, plus factors of weight `m`} between
          colours) whose restriction `I^- -> I^+` is injective, with a unit
          minor `(-1)^m I`.  Hence `Phi = tr` of an orthogonal projection on
          `I^+`.  This is unbounded in labels and factor count, inside the
          known positive region (part `>= sum - 1`).
        - Proposition 2: every equivariant endomorphism supported on at most
          two factors of a pair-free list commutes with `tau`.
        - Theorem 2: for `(-n)^2 (+1)^(2n)`, the smallest support of a
          `tau`-odd equivariant endomorphism is `ceil(n/2)`.  So no fixed
          support bound works uniformly (knob: locality); an injection must
          move blocks that grow with the labels.
        - Obstructions to extending unchanged (knobs: maximal-spin channels
          only; pure-colour targets; a single stage of singlet transfers).
        - Open (first unresolved formula): for general two-minus lists,
          all invariant channels and mixed-colour targets, `det((J(t)|I^-)^*
          (J(t)|I^-)) != 0` for some weights `t`.
      - FM-STR5 (luna_max_mars; `fm39/str5_representation_tests.py`, rerun: same
        output).  Exact tests for the representation routes.
        - Invariant bases are explicit as Clebsch-Gordan path pairs.  The
          plain factor swap has equal eigenspaces; the relevant involution is
          the sign-twisted one.
        - Three canonical singlet transports have rank deficits, or vanish
          on `(1^4, 2^2)`.
        - Groups: in `Sp(4)`, `h_n` is genuine but `S_p = V_(p,0) - V_(p-1,1) +
          V_(p-2,0)` (`p >= 2`).  `SO(5)`, `Sp(2r+2)`, `SO(2r+3)`, `SU(3)`,
          `G_2` and `T^2` fail as Haar models for `r >= 2`.
        - Ginibre torus: plus factors are positive definite; minus factors
          have coefficient -1 at `(n,-n)`.
      - FM-SEC177 (luna_max_venus; `fm39/sec177_42_factor_noflip_check.py`;
        main-agent rerun: all 861 flips negative, TopPair monotone;
        independent check of the 18-factor member with
        `sec166_flip_single.cpp`).  NO-FLIP LISTS HAVE UNBOUNDED FACTOR COUNT
        (in all tested cases).
        - The all-minus consecutive runs `B = (-1, ..., -k)` with `sigma p` near
          `+-(k+1)` have no flip descent for every tested `k` from 17 to 41.
          At `k = 41`, `sigma p = -43`: 42 factors, `W = 861`, `delta = 409`,
          largest flip `-1.39e43`.
        - TopPair is monotone on every one; the child/parent ratio is
          `<= 0.0072` (the ratio rises slowly with `k`).
        - So a factor-count closure, plus a many-factor flip theorem,
          cannot cover the cone: hard lists exist at every factor count.
          This confirms the reallocation to structural mechanisms.  Exact
          class screens of doubled, tripled and quadrupled runs found no
          no-flip pattern.
      - FM-CHK108 (luna_max_eris, fresh code;
        `fm39/chk108_ten_factor_check_part1.py` .. `part4.py`, all rerun:
        PASS).  ACCEPT FM-MECH169 (ten factors) on all four items:
        - the reductions;
        - the profiles `H^(9)`, `H^(10)` (shape counts 8,930 and 6,336
          reproduced);
        - the payment tables (maxima below 15/16);
        - the finite box, checked against its cutoff with a one-million-
          signing exact sample (no negatives) and six boundary profiles, and
          all 1,025,024 signings of ten-factor words with labels 1..5.
        The full box itself was rerun by the main agent.
      - FM-STR8 (astra_max_vulcan; `fm39/str8_block_transfer_complex.py`, rerun:
        PASS).  A complex with Euler characteristic `Phi`; homology not yet
        concentrated.
        - Obstruction 1 (knob: the chain groups `C_k` graded by `|S cap M|` with
          adjacent-degree maps only): at `(-1,-1,-2,-3,-3,-4)`, `C = (11,0,4,10,
          4,0,11)` admits no even concentration, even after a parity shift.
        - Construction (proved): order the minus positions.  For each block
          `R` with `|R cap M|` odd, let `Q_R` project onto `Inv(V_R)` and move the
          singlet block to the other colour, with exterior-contraction signs.
          `d = sum_R c_R Q_R` is odd and `d^2 = 0` for any scalars `c_R`, so `Phi =
          sum_k (-1)^k dim H_k`.
        - Homology (with `c_R = 1`): `(1,0,1)` for `(-1,-1,+1,+1)`; `(6,0,4,0,4,0,6)`
          at the word above; `(10,0,45,10,45,0,10)` for `(-2)^6`, with an explicit
          surviving odd cycle.
        - Theorem (knob: nested-cut singlet transfers only): on `(-2)^6`, every
          square-zero odd differential built from scalar multiples of
          nested-cut singlet transfers, in either direction, has `rank O +
          rank I <= 10 < 20` (the triple-cut images span only 5 dimensions).
        - Together with FM-STR1's single-stage obstruction: singlet transfers
          alone cannot give an injection or a concentrated complex.  The next
          ingredient is transfers through nonzero intermediate fusion
          channels (recoupling with 6j coefficients) or exchanges across a
          cut.
      - FM-STR4 (luna_max_pluto; `fm39/str4_channel_identity_mining.py`, rerun:
        ALL CHECKS PASS).  First exact positive descent identities, found by
        exact LP.
        - For `B = (1,2,5,9,11,12)`, `p = 14` (TopPair `(9,11)`) and for `B =
          (1,3,4,6,7,9)`, `p = 12` (TopPair `(7,9)`):
          - the TopPair difference equals a NONNEGATIVE rational combination
            of shorter channel-fusion children `Phi_((a,b)->c)` (two factors
            replaced by one channel `c`), as an identity for all 64
            signings;
          - e.g. the first: `Delta_TP = (212/1331) Phi_((1,2)->3) + (271/1331)
            Phi_((1,11)->10) + ...`, 7 terms, L1-optimal over 149 child
            columns;
          - so TopPair descends there, by induction, without any flip
            hypothesis.
        - Such an all-signing identity cannot exist where TopPair fails for
          some signing (the 126 census failures): flip terms must enter.
        - Knob: an aggregate-`Q` coefficient system is infeasible (Farkas
          certificate), while the channel-resolved system succeeds.
      - FM-STR7 (astra_max_minerva; `fm39/str7_folded_involution.py`, rerun:
        FM-STR7 PASS).  A sign-reversing involution on the diagram model.
        - Model: signed diagram pairs = fusion paths in pairs of coordinates,
          one Clebsch-Gordan step per factor; the noncrossing diagrams are
          recovered by stack matching.
        - Folded involution (proved): if every signed label has even
          multiplicity, order the word as `H H`.  At each midpoint pair
          half-paths of opposite sign lexicographically and reflect.  This
          is sign-reversing with only positive fixed points, so `Phi >= 0`.
          It covers every fundamental-only word, at every length.  It
          changes the uncolored matching, so it is not the failed overlay
          argument.
        - The same construction realizes flip descent combinatorially.
        - Obstructions:
          - Preserving the intermediate channel fails on an actual no-flip
            residual: every cut with `>= 2` factors per side has a negative
            fiber.
          - Recoloring a bounded number of factors cannot work (Hall: demand
            8 > capacity 4, at every odd scaling); a uniform family forces
            unboundedly many factor reassignments.
          - Knobs: channel preservation; bounded recoloring.
        - Open: for some pair, an injection from the `-2 D_uv` negative mixed
          remnants into the `2 sum_c A_c0(C)` positive pure remnants, using
          `Phi = 2(sum_c A_c0(C) + D_uv)`.  It must change the channel and may
          move arbitrarily many plus factors.
      - FM-STR6 (luna_max_neptune; `fm39/str6_noflip_census_laws.py`; main-agent
        rerun: same counts).  Census laws for the 75,532 no-flip lists at `W <=
        52`:
        - label multiplicity `<= 4`; at most two labels 1; the number of odd
          factors is never 2 (that case is two-odd fusion);
        - the minimizing signs are not determined by the labels (e.g. eight
          no-flip signings of `(1..7,9,10,11)`);
        - sign sequences are neither necessarily monotone nor alternating.
        - None of these laws isolates no-flip lists: 120,743,304 residual
          lists satisfy all of them.
        - Its conjecture `D <= 10`, `L <= 16` (`D` distinct labels) is REFUTED by
          FM-SEC177: `(-1, ..., -41)` with `sigma p = -43` is no-flip with
          `D = 42` and `L = 42`.  The census bound is an artifact of `W <= 52`.
      - FM-STR1b (astra_max_ceres; `fm39/str1b_two_stage_injection.py`, rerun:
        FM-STR1b PASS).  TWO-STAGE compositions of transfers.
        - Proposition 1 (kill, knob: one stage of all-one-colour singlet
          transfers): for `(-1)^2 (+3)^6`, every one-stage `J(t)`, with any
          weights and all channels, has rank `<= 128` on the 160-dimensional
          `I^-`.  An exchange-symmetry count (skew source 40 > skew target 34)
          shows this before any computation.  FM3 holds there (`Phi = 786`).
        - Proposition 2: a composition (block transfer, then a pair singlet
          transfer, with endpoint weights) is injective on that example (an
          80 x 80 minor nonzero mod the prime 65521).
        - Theorem 3: for four minus labels 2 and a plus background with
          `m_2(S) = 0` for every proper nonempty plus subset (condition PB),
          the explicit two-stage operator `sum_(i<j) K_ij F_i` is injective on
          `I^-` (`dim I^- = 8b`), via `SO(3)` identities on `V_2 = C^3` (a 4 x 4
          minor of determinant 4; each `Psi_ij` an isometry).
          - This covers every `(-2)^4 (+q)(+t)` with `q, t != 2`, and
            `(-2)^4 (4, 8, ..., 2^L, 2^(L+1)-2)` for every `L`: unbounded
            factor count.
          - So `Phi = tr` of a positive projection there.
        - Proposition 4: compositions also make `(-2)^6` injective (rank
          20/20), escaping the FM-STR8 nested-singlet bound.  The second
          transfer crosses the original colour cut.
        - Convergence: injection (FM-STR1/1b), complex (FM-STR8) and the
          identity-depth result (FM-STR3) all need TWO levels of recoupling.
          Open: a general multi-stage rule (e.g. along a fusion tree) for
          two-minus and `2r`-minus lists without condition PB.
      - FM-STR6b (luna_max_neptune; `fm39/str6b_mechanism_survey_checks.py`,
        rerun: same output).  Survey of known mechanisms for Ginibre /
        GKS-II-type inequalities, translated to FM3.  Citations are from the
        agent's memory, with DOIs as given; not checked against the
        literature by the main agent.
        - Mechanisms covered:
          - Ginibre (CMP 16, 1970) duplicate variables / positive
            definiteness;
          - Sylvester (CMP 73, 1980), a counterexample to the
            positive-definiteness inequality for `O(N)`, `N >= 3`, which is
            not a counterexample to Griffiths inequalities or FM3;
          - Dunlop-Newman (CMP 44, 1975); Dunlop (CMP 49, 1976);
            Kunz-Pfister-Vuillermot (1975/76);
          - Messager-Miracle-Sole-Pfister (CMP 58, 1978); Bricmont (JSP 17,
            1977); Wells (bibliographic source UNCERTAIN);
          - Newman (CMP 41, 1975), Lee-Yang;
          - Simon-Griffiths, Aizenman, Brydges-Frohlich-Spencer/Sokal
            (random currents and walks);
          - Frohlich-Israel-Lieb-Simon (CMP 62, 1978), reflection
            positivity.
        - Exact checks: all 6,469 signed multisets with labels 1..4 and
          length `<= 8` are nonnegative.  Factorwise positive definiteness
          fails already at `(-1,-1)` (coefficient -2).
        - Recommendation: (1) global SU(2) spin-network switching (current
          switching with full fusion trees, nonzero channels and
          multi-generation recoupling, consistent with FM-STR1b/3/7/8); (2) a
          whole-word Gram or reflection-positive cut construction.
      - FM-STR8b (astra_max_vulcan; `fm39/str8b_channel_complex.py`, rerun: PASS,
        481 profiles).  A complex using ALL intermediate channels; homology
        concentrated on the whole test box.
        - Setup: sort the factors by signed class and split alternating
          positions into half-words `A, B`.  On each half's super-
          representation, match odd to even copies of the same `SU(2) x SU(2)`
          type, giving `d_A`, `d_B`.  Put `d_0 = d_A (x) 1 + Gamma_A (x) d_B`.  Then
          `d_0^2 = 0`, and with `f_A(alpha) = e_A(alpha) - o_A(alpha)` (signed
          multiplicity of type `alpha`):
            `dim H_even = sum_alpha max(f_A f_B, 0)`,
            `dim H_odd = sum_alpha max(-f_A f_B, 0)`,
          so `Phi = sum_alpha f_A(alpha) f_B(alpha)`.
        - Proved: concentration whenever `f_A, f_B` have compatible signs in
          every channel, in particular for every doubled word `Gamma Gamma`
          (`dim H_even = sum f^2`), uniformly in labels and factor count.
        - On `(-2)^6`, the nonzero channel `(2,2)` supplies 18 of 20
          cancellations, giving homology `(100,0)`.  This escapes the
          FM-STR8 nested-singlet bound.
        - Exchange correction: `delta = P J P`, with `J_TS = [prod_(i in S cap T)
          (i+2)] i_T^* i_S` from odd to even cuts and `P` the projection onto
          the unmatched channel spaces, fixed before testing.  Then `d = d_0 +
          delta` has `d^2 = 0` and `dim H_odd(d) = dim H_odd(d_0) - rank delta`.
        - Knob: channel preservation, e.g. `(-1)^3 (-2)^3 (+3)` leaves 4 odd
          dimensions for any channel-preserving differential; the
          correction removes them.
        - Exact screen: all 481 pair-free signed multisets with labels `<= 3`,
          `<= 8` factors, even minus count.  `d_0` already concentrates in 471;
          the 10 others are fully corrected by `delta`.  `H_odd(d) = 0` in
          every case, so `Phi = dim H_even` there.
        - Open (FM3 needs only this): `delta` restricted to `H_odd(d_0) ->
          H_even(d_0)` is injective for every pair-free residual word.
      - FM-CHK110 (luna_max_neptune, fresh code; `fm39/chk110_structural_theorems_check.py`,
        rerun: PASS).
        - ACCEPT FM-STR1b Theorem 3 under PB, with both corollary families:
          66 label pairs, `A_L` for `L = 2..8`, rank `8b`, the minor of
          determinant 4, the Schouten and `delta`-contraction identities.
        - ACCEPT FM-STR7 Proposition 2 (the folded involution).  Each
          endpoint contributes `(P_v - N_v)^2` fixed pairs, so `Phi = sum_v (P_v
          - N_v)^2`.  Checked on all 495 even-multiplicity words with labels
          `<= 4` and length `<= 8`.
      - FM-STR3b (astra_max_juno; `fm39/str3b_fusion_depth_lp.py`, rerun: ALL
        EXACT CHECKS PASS).  Fusion depth by LP duality.
        - Setup: product-variable LP `Y(A,B) ~ m(A)m(B)` closed under fusion to
          depth `k` (a linear relaxation of the quadratic problem).
        - With only unit/singleton invariants fixed, the tested depths are 1
          or 2.  Fixing the pair values `m(a,b) = delta_ab` gives ONE generation
          on every tested root: a six-factor root, a seven-factor root, and
          all four strict no-flip signings of labels 1..8.
        - Every extracted dual certificate has ZERO flip coefficients: it is
          an unconditional positive TopPair identity there (e.g. three
          shorter forms for the 1..8 signing A).
        - Depth zero always fails: setting only the root invariant to zero
          makes TopPair negative.
        - The FM-STR3 one-generation obstruction is explained: it violates
          the derivative of the pair normalization `m(n,n) = 1`.  So it does
          not show a need for deep 6j recoupling.
        - Undetermined: the (F1) root and the eleven-factor root (complete
          recoupling bounds the depth by 14 and 10).
      - FM-STR4b (luna_max_pluto; `fm39/str4b_depth_one_cone_tests.py`, rerun:
        same output).
        - Exact Farkas kill (knob: all-signing identities in the depth-one
          channel dictionary without recoupling relations, even with any
          flip terms): `B = (1,2,3^3,4^2,5^3)`, `p = 15` has TopPair positive on
          all 32 signings but no such identity.  A dimension-, 6j- or
          symmetric-weight normalization in that dictionary is therefore
          impossible at all-signing strength.
        - `B = (1^18, 2, 3^4, 4)`, `p = 6`: the minimum flip support is exactly
          1 (`D_(1,1)`).
        - Main-agent caveat: its "no-flip-only" identity on the first profile
          (one child, two signings) is fitted pointwise and carries no
          structural information.  Juno's product-variable certificates,
          valid for every table satisfying the constraints, are the
          meaningful kind.
        - The `W <= 36` census map was not completed.
      - FM-STR7b (astra_max_minerva; `fm39/str7b_height_transport.py`, rerun:
        FM-STR7b PASS).  Folding with a remainder.
        - For `Lambda = H H R`: `Phi = sum_(u,v) r_H(u,v) r_(HR)(u,v)`, where `r` is the
          signed half-path count, i.e. the two-colour table of the half-word.
        - Proposition 2: after midpoint and same-height cancellation, a
          channel-changing involution with only positive fixed points exists
          IF AND ONLY IF the height prefixes `P_T = sum_(u+v <= T) r_H r_(HR) >= 0`
          for every `T` (match each negative remnant to a positive one of
          lower height).  This changes channels, heights and whole
          half-paths.
        - Exact obstructions (knobs: midpoint preservation, one-remainder-
          per-side, height preservation, pooling two adjacent heights,
          PSD on unrestricted half-vectors, "two remainders force a flip").
          E.g. `H = (-6,-8,-10)`, `R = (-2,-4)` has `Phi = 8,226` and all 28 flips
          negative.
        - Induction on `|R|` needs repair: `|R| = 1` occurs (one unpaired even
          plus label); restoring pair-freeness after fusion can raise `|R|`
          from 2 to `k+1`.  Repair: allow all signed-class words as children
          and postpone pair reduction.
        - Census: every signed multiset with labels 1..4, `<= 10` factors, `|R| =
          2` or 4 (1,860 + 280 pair-free; 5,940 + 6,270 all classes) has every
          prefix `P_T >= 0`, and so do 6,000 seeded stress cases.
        - Open: the prefix positivity `sum_(u+v <= T) r_H r_(H, eps a, eps b) >= 0`
          for all `T` (its `T = infinity` case is FM3 for two remainders).
      - FM-STR1c (astra_max_ceres; `fm39/str1c_two_cut_injection.py`, rerun:
        FM-STR1c PASS).
        - Theorem 1 (structural no-go, knob: naturality under permutations
          of equal signed factors): for a two-minus list with six plus
          factors, no `S_2 x S_6`-equivariant injection `I^- -> I^+` exists, with
          any channels, colourings or number of stages.  Exact character
          calculation: a nine-dimensional `S_6` component is killed by every
          equivariant map.  So any injection must break symmetry (an
          ordering or a distinguished factor), as the lexicographic folding
          of FM-STR7 does.
        - Theorem 2 (proved): under a two-complementary-cut support
          condition (TC), `T_A + T_B` is injective on `I^-`, with an exact energy
          identity: an isometry if `n != m`, and if `n = m = q`, `det(J^*J) =
          [q(q+2)/(q+1)^2]^(2b)`.  Corollary: an explicit family with `L + 5`
          factors for every `L`.
        - Obstructions: fixed-prefix first stages (rank `<= 2`); corrections made
          only of plus-pair singlet transfers can annihilate whole source
          channels.  Keeping a direct branch repairs the tested cases
          (`3^6, 3^5 5, 3^4 5^2, 3^3 5^3`: rank 160/160, minors nonzero mod
          65521).
        - Open: combine the cut maps `U_S` (with `sum_S U_S^* U_S = I`) inside the
          actual `I^+` for every two-minus list.
      - FM-CHK109 (luna_max_eris, fresh code; `fm39/chk109_eleven_factor_partial_check.py`,
        rerun: PASS).  FM-MECH173 (eleven factors), partial acceptance.
        - Accepted:
          - the payment identity and reductions;
          - the profile catalogs (57,920,058 normal and 27,446,154 all-even
            shapes);
          - the finite-box cutoff;
          - exact samples: 512 envelope profiles and boundary profiles, plus
            all 745,472 signings of eleven-factor words with labels 1..4.
        - REPAIR NEEDED (knob: seven-small-prefix budget coverage): the
          universal charge bounds for the 31,802 claimed-good seven-label
          prefixes were not independently regenerated.  This is a gap in
          the check, not a counterexample.
        - The strata row for eleven factors stays "independent check
          pending".
      - FM-STR7c (astra_max_minerva; `fm39/str7c_fundamental_prefix_positivity.py`,
        rerun: FM-STR7c PASS).  Height-prefix positivity, proved for
        fundamental half-words.
        - Lemma 2: the height-truncated fundamental row `Pi_T (x+y)^h` is the
          restricted character of a GENUINE `Sp(4)` module `C_(h,T)`.  The E-basis
          coefficients are nonnegative; the boundary layer `u+v = tau+1` needs
          an explicit ballot inequality (it is not a highest-weight
          truncation of the tensor power).
        - Theorem 3: for `H = (+-1)^h` and every admissible remainder `(eps a, eps
          b)`, `P_T >= 0` for every `T`.  With `eps = -1` the certificate is an
          `Sp(4)` pairing of `C_(h,T)` with `V^(x)h (x) Sym^(a-1)V (x) Sym^(b-1)V`.
          So the FM-STR7b involution is fully proved for all words
          `(+-1)^(2h) (eps a)(eps b)`.  Their scalar FM3 was already known; the
          new content is the mechanism.
        - Proposition 4: the height cutoff is not the finite-level `SU(2)_k`
          corner.  The exact relation keeps a projector; neither equality
          nor a one-sided comparison holds (exact tables).  Knobs: equality,
          one-sided comparison, commuting compressed generators.
        - Proposition 5 (knob: positivity at each `Sp(4)` highest weight
          separately): for general `H` the termwise certificate fails, e.g. `H
          = (+2)`, `R = (-3,-5)`; cancellation between highest weights is
          needed.  306,036 prefixes over 25,560 one-general-label inputs:
          no negative prefix.
        - Open: `H = (+1)^h (eta n)`, one general label.
      - FM-CHK111 (luna_max_neptune, fresh code; `fm39/chk111_str8b_independent_check.py`,
        rerun: PASS, about 4 min).  Independent exact-rational check of the
        FM-STR8b all-channel complex.  ACCEPT on the stated scopes.
        - `d_0^2 = 0` (tensor differential identity) and the homology formula
          for `d_0` are reproduced.  With the correction `delta = PJP`,
          `delta^2 = 0` and `d_0 delta = delta d_0 = 0`, so `d^2 = 0`.
        - Screen of 481 profiles: 88 have nonzero odd chains and 10 have
          nonzero `H_odd(d_0)`.  All 10 correction matrices have full
          column rank, so `H_odd(d) = 0` on the whole screen.  Largest
          case: `+1,-2^2,+3^5`, `H_0 = (418,10)`, rank 10.
        - Random profiles of length 9 (one with `H_0 = (638,2)`, rank 2) and
          12 seeded samples each at lengths 9 and 10: all pass.
        - Agrees with the FM-STR8c stress log (4,521 pair-free profiles,
          labels `<= 4`, length `<= 10`; prime-field full-rank minors).
        - Open (knob: uniformity over labels and factor count):
          `rank(PJP : H_odd(d_0) -> H_even(d_0)) = dim H_odd(d_0)` for every
          pair-free list.
      - FM-STR9 (astra_max_ceres; `fm39/str9_two_cut_injectivity.py`, rerun:
        FM-STR9 PASS, 8 s).  Injectivity of the FM-STR8b correction
        `delta = PJP` on a structural two-minus class.
        - Proposition 1: for two-minus lists the unmatched odd channels are
          given exactly by fusion operators `L_n, R_n` acting on the
          unsigned two-colour table `G_C`.  One minus per half: both tables
          antisymmetric, odd channels in reflected pairs.  Both minuses in
          one half: odd channels are where `G_D > 0` and `f_A < 0`.
        - Theorem 2: let `Lambda = (-n)(-m) prod (+a_i)`.  Assume (TC): only
          one complementary pair of plus-block cuts couples, i.e.
          `b_n(R) b_m(I - R) = 0` unless `R = P_0` or `R = Q_0`.  Assume (NS):
          no nonempty sublist inside one half has an invariant.  Then
          `delta` is injective on `H_odd(d_0)`, with the energy bound
          `|PJP xi|^2 >= c |xi|^2`, where
          `c = min_j (a_j - b_j)^2 / (2 + a_j^2 + b_j^2)` (an extra factor
          `q/(q+1)` when `n = m = q`), and an explicit Gram determinant.
          `a_j = b_j` is excluded by Bertrand's postulate.  The
          certificate uses only the pure-colour targets `(t,0)` and `(0,t)`.
        - Corollary 3: an explicit family `(-1)^2 (+u)(+v)(+w) H_l` with
          `l + 5` factors for every `l`, and `dim H_odd(d_0) >= 2 r_0 (l - 1)
          > 0`.  The correction is injective there, so FM3 holds on the
          family by a positive-operator trace.  First uniform class with
          unbounded factor count settled by the complex.
        - Proposition 4 (knob: correction image restricted to pure-colour
          targets): for `(+2)^3(+3)^2(-4)(+5)^4(+6)^2(+7)(-9)`,
          `H(d_0) = (56,676,672, 17,240,634)`.  The pure targets have total
          dimension 15,925,636, leaving a kernel of dimension at least
          1,314,998.  This does not bound the full mixed-colour correction.
        - Open: `<xi, J^* Q J xi> > 0` for every nonzero `xi` in `H_odd(d_0)`
          with `R J xi = 0`, where `R` projects to the pure-colour part of
          `H_even(d_0)` and `Q = P_even - R`.
      - FM-STR3c (astra_max_juno; `fm39/str3c_fused_child_obstruction.py`,
        rerun: ALL EXACT CHECKS PASS, 41 s).  Pair-normalized depth-one
        identities for the TopPair difference.
        - The four no-flip signings of labels 1..8 (negative sets `{1,5}`,
          `{3,7}`, `{2,4,6,8}`, all eight) have exact identities
          `Delta_TP = q + sum gamma m(A)m(B) + sum lambda R(A,B;i,j)` with
          `gamma >= 0` and depth-one fusion relations `R`.  These need no
          shorter FM3 values and no flip terms.  Example (all eight
          negative): `Delta_TP = 472`, `q = 155`, 545 fusion terms, 155
          product terms.  The coefficients come from rational dual
          extraction; no closed formula yet.
        - Depth ranges: eleven-factor root `Delta_TP = 191,001`, `1 <= k* <=
          9`; F1 `Delta_TP = 3,914,630`, `1 <= k* <= 13`.  Depth zero fails
          for both.
        - Obstruction (knob: the shorter-form dictionary; pair
          normalization and all shared-background depth-one recouplings
          kept).  For F1 there is a nonnegative integer direction `v`
          (1,090 nonzero entries, zero on unit, singletons and pairs) that
          satisfies all 10,448 depth-one scalar fusion equations.  It has
          `d_v D_uv <= 0` for every root flip, and `d_v Phi >= 0` on every
          proper deletion child and every fused word with at most ten
          factors (700,367 signed profiles, exact).  But `d_v Delta_TP =
          -297,813`.  So any such identity for F1 must use a genuinely
          fused child with 11 to 14 factors.  Example of an excluded
          longer form: `(+1,+1,+2,-3^3,-4^3,-5^4,+11)` (channel
          `(3,8) -> 11`), with `d_v Phi = -5,069,906`.
        - What it gives the cone: no new stratum.  The induction behind
          TopPair must consume positivity of long fused children, not just
          products of invariant tables.
        - Open: on F1, `Delta_TP + sum alpha D_uv = sum beta Phi(Q,eta) + sum
          gamma m(A)m(B) + sum rho R_J` with `alpha, beta, gamma >= 0`, where
          some `beta` is positive on a fused word with `11 <= |Q| <= 14`.
      - FM-CHK109b (luna_max_eris, fresh code; `fm39/chk109b_prefix_budget_check.py`,
        rerun: CHARGE PREFIX SCREEN PASS, 35 s with 24 threads).  ACCEPT
        FM-MECH173 item 3, which closes the gap left by FM-CHK109.
        - All 31,824 sorted seven-label prefixes were regenerated (8,191,650
          feasible sign/category cases in mode 0, 27,090 in the all-even
          mode).  31,802 satisfy the exact `2^40` charge bound.  Worst:
          `(1,1,1,2,2,3,4)`, numerator 1,093,158,045,286, margin
          6,353,582,490.
        - The 22 exceptions and their exact maxima match the claimed list.
          Each goes to the finite phase with cutoff `q <= 2D <= 44`.
        - The eleven-factor strata row is now marked independently checked:
          FM-CHK109 covers the reductions, profiles, cutoff and samples, and
          FM-CHK109b covers the prefix budgets.  The finite box itself is the
          main-agent rerun of the FM-MECH173 reproducer.
      - FM-STR7d (astra_max_minerva; `fm39/str7d_sp4_truncation.py`, rerun:
        FM-STR7d PASS, 27 s).  The `Sp(4)` height-prefix mechanism beyond
        fundamentals.  Not closed for an arbitrary general label.
        - Lemma 1: the height cutoff `Pi_T` sends every genuine `Sp(4)`
          character `chi_(alpha,beta)` to a genuine one: zero if `T <
          alpha - beta`, else `chi_(l+J,J)` with `l = alpha - beta`, `J =
          min(beta, floor((T-l)/2))`.  For antisymmetric rows, `Pi_T(dB) = d
          J_(T-1) B`.  This builds a new module; it is not a submodule
          filtration of the original (`Pi_0(V (x) V) = 2 . 1`).
        - Theorem 2: if the half-word row `F_H` lifts to a GENUINE `Sp(4)`
          character `A`, then every prefix is nonnegative for either common
          remainder sign: `P_T(H; -a, -b) = 2 <L_T A, A K_a K_b> >= 0`.  The
          same holds for remainders `d^2 B` with `B` genuine.
        - Corollary 3: prefix positivity for `(+-1)^h (+2)^r` with `h >= r`
          (from `s S_2 = chi_(3,0) + chi_(1,0)`).  This gives the FM-STR7b
          matching for `(+-1)^(2h) (+2)^(2r) (eps a)(eps b)`.  The scalar
          values were already known.
        - Proposition 4 (three stronger statements fail exactly):
          - knob: a fixed number of fundamentals makes `S_n` genuine.
            `[chi_(n-1,h+1)] s^h S_n = -1` for `n >= h + 2`.
          - knob: arbitrary genuine `Sp(4)` multiplicities in place of the
            actual tensor power.  `G = chi_(5,0) + 3 chi_(4,3)`, `F = G
            S_5`, `R = D_2 D_4` gives `P_0 = -4`.
          - knob: certify each constituent separately.  `chi_(3,2) D_3`
            with `R = D_2 D_8` has `L_4 = -2`, while the actual word
            `s^5 D_3` has `P_4 = 32,180`.  Cancellation between
            constituents is needed.
        - Proposition 5: an explicit quartic multiplier,
          `[d^2 s^t : chi_(alpha,beta)] = m_(t+2)(alpha,beta) Q_t(A,B) / D_t`,
          `Q_t = A^2 B^2 - (3t+13)(A^2+B^2) + 15t^2 + 120t + 241`.  It is
          negative at the balanced top weight for every `t`, so `d^2 s^t` is
          never genuine.
        - Proposition 6 (knob: the fusion correction itself nonnegative):
          the induction on the remainder size `|R|` still needs a signed
          estimate.  `Phi(HH,-2,-4) = 30 - 10` for `H = (+1,+3)`.
        - Screens: 408,776 prefixes on 4,955 one- or two-general-label
          inputs, and 148,898 more with `h <= 200`.  No negative prefix.
        - Open: for `H = (+1)^h (-2)`,
          `<J_(T-1) V^(x)(h+1), V^(x)(h+1) d^2 K_a K_b>_Sp(4) >= 0` for every `T`
          and every `2 <= a < b` with `a = b mod 2`.
      - FM-STR8c (astra_max_vulcan; `fm39/str8c_stress_and_tripled_classes.py`,
        rerun: PASS completed verifier, about 3 min).  Stress test of the fixed FM-STR8b
        correction, plus a uniform proof for two tripled classes.
        - Box: all 4,521 pair-free profiles with labels `<= 4`, at most ten
          factors, even minus count.  3,990 have `O = 0` already; the other
          531 all have `rank delta = O` (269 matrices after odd reflection;
          largest rank 4,918).  Ranks are certified by full-rank minors
          modulo 1,000,003 with every denominator checked.
        - Named lists, `H(d_0) -> H(d)`, all concentrated:
          - `(-1,+2,+3,+4,-5,+6,+7,+8)`: `(1,172, 216) -> (956, 0)`;
          - runs `(-1,...,-k, (-1)^k p)`: k = 5 `(40,4)`, k = 6 `(196,12)`,
            k = 7 `(1,080,100)`, k = 8 `(7,256,822)`, k = 9
            `(53,872, 7,664) -> (46,208, 0)`.
        - Not completed (dimensions only): `(+1)^2(-2)(+3)^4(-4)^3(+5)^4(+8)`
          with `(E,O) = (8,437,374, 286,632)`; the even run `(-40,+42,-44,+46,
          +48,...,+62)` with `O = 2,445,554,811,680`.
        - Proposition 2 (proved, all labels): for `(-a)^3(-b)^3`, `H_odd(d) =
          0`.  `O = 2` exactly when `a, b` are even and `a < b <= 2a`, else
          `O = 0`.  In the `O = 2` case the retained minor is diagonal with
          `|det| = 6 U(s,t)^2 > 0` (Racah sum with a single index).  The
          exchanges go across nonnested cuts through channel `a+b`.  These
          words were already inside known FM3 strata; the new content is
          the mechanism.
        - Observation: the odd fraction `O/E` grows along the runs (10%, 6%,
          9%, 11%, 14% for k = 5..9).  It is about 30% on the FM-STR9
          fourteen-factor list.
        - Open (knob: simultaneous control of several incompatible
          channels): `ker(P_even J P_odd) = 0` for every pair-free
          even-minus list.
      - FM-STR4c (luna_max_pluto; `fm39/str4c_depth_lp_coverage.py`, rerun:
        same summary, about 10 min; reads the gzipped `W <= 40` census logs
        `fm39/sec166_census_w40_noflip.log.gz` and `..._tpfail.log.gz`).
        Coverage map of the pair-normalized
        depth LP on the no-flip rows.  Finite screen, not a uniform result.
        - The 5,430 no-flip rows group into 2,014 label-multiset domains.
          Screened: all 153 six-factor domains and the first 323
          seven-factor domains (1,476 signings).
        - Depth one certifies 425 of 476 domains (89.3%).  Each of the other
          51 has an exact depth-one separator (46 negative points, 5 negative
          rays) and is certified at depth two.  No certificate needs flip
          multipliers or shorter-list multipliers at depth two.  Example
          `(2,3,5,7,8,9)`: depth-one optimum `-1`, actual drop 61, depth-two
          bound 61.
        - The 126 TopPair-failure rows (63 multisets, 976 pair-free signings)
          have no no-flip signing, so the flip route applies to all of them.
        - Untested: 1,538 domains (195 seven-factor, all with `>= 8` factors).
        - Caution (main agent): a depth-two bound equal to the actual value
          may mean the depth-two relations pin every invariant that enters.
          Then the certificate is only an evaluation.  FM-STR4d checks this
          first.
      - FM-STR9b (astra_max_ceres; `fm39/str9b_mixed_colour_targets.py`,
        rerun: PASS, about 15 s).  Mixed-colour targets for the FM-STR8b
        correction.
        - Proposition 1 (highest-weight slice Gram identity): projecting a
          large label `V_n` to its highest-weight coefficient scales every
          invariant Gram entry by the same factor `n+1`.  Gram
          computations with one large label then live in a fixed finite
          coordinate space, independent of `n`.
        - Theorem 2 (uniform in `n`): for `Lambda_n = (+1)^5(-2)(+n)(-n-1)`,
          `n >= 4`, `H(d_0) = (162, 22)`.  The pure-colour map has a
          two-dimensional kernel; two retained mixed targets (half-channel
          `(0, n-1)`) repair it, with `<xi, J^*QJ xi> >= |xi|^2` on `ker RJ`.
          All Gram entries are identities in `Q(n)`.  Corollary 3: the same
          for the odd reflections (six and eight minus factors).  These
          words were already inside known strata.
        - Proposition 4: `(+1)^7(-2)(+3)(-4)` has `H(d_0) = (1,678, 336)`, pure
          rank 150, full rank 336.  Mixed targets repair a 186-dimensional
          kernel (modulo 65521, denominators checked).
        - Theorem 5 (target capacity; knob: a fixed number of mixed
          targets).  On `Lambda_t = (+1)^(2t)(-2)(+4)(-6)` the odd fraction
          `h_o/h_e` tends to a limit `rho` with `1/163 < rho < 1`, both of
          order `16^t / t^5`.  The pure targets have dimension at most `210 .
          4^t`, so `dim ker(RJ) / h_o -> 1`.  Method: the half-word tables
          have Gaussian scaling limits, `f_A ~ Psi h`, `f_B ~ h`, with `h =
          xy e^(-(x^2+y^2)/2)` and `Psi = (x^2-y^2)^2 - 4(x^2+y^2) + 12`.  The
          limiting Euler characteristic is `int x^2 y^2 e^(-r^2) Psi = 3 pi /
          16 > 0`.  Exact ratios: 0.06-0.07 (t = 4) up to 0.23-0.24 (t = 32).
        - So the correction must handle mixed target spaces of growing
          dimension, almost all of the odd homology on this family.
        - Open: for `Lambda_t`, uniformly in `t >= 5`, `<xi, J_t^* Q_t J_t xi>
          > 0` for every nonzero `xi` in `H_odd(d_0)` with `R_t J_t xi = 0`.
      - FM-STR3d (astra_max_juno; `fm39/str3d_complete_block_obstruction.py`,
        rerun: ALL EXACT CHECKS PASS, 9 s).  Complete fusion blocks are
        insufficient for TopPair on F1.
        - Lemma 1: `Q_(uv,eta) = (Phi(L,eta) + Phi(L,eta^(uv)))/2 = sum_c
          Phi(L - {u,v}, eta_u eta_v c)`: a complete block has equal
          weights over its CG channels and the inherited sign.
        - Proposition 2 (knob: complete channel sums instead of individually
          weighted channels).  For F1 = `(+1)^2(-2)(+3)^4(-4)^3(+5)^4(+8)`
          there is an integer direction `v` (3,036 nonzero entries).  It
          fixes all invariants with at most three factors and all
          odd-weight invariants, and satisfies all 10,448 depth-one
          scalar fusion equations.  Along `v`, `Delta_TP` has derivative
          -77,799,242,809, every `4 D_uv` is negative, and every re-signed
          complete block `2 Q_(uv,eta)` is `>= 0` (all 19,592 sign-count block
          pairs).  So no depth-one identity with nonnegative complete-block
          coefficients, any re-signings and any root flips exists.
        - Proposition 3: an individual channel detects `v`.  The p-incident
          top channel `(4,8) -> 12`, child `(+1)^2(-2)(+3)^4(+4)^2(+5)^4(-12)`,
          has derivative -2,203,671,549,030, while its complete block has
          +1,599,180,810,932.
        - Proposition 4: fusion induction is well founded in the order
          (total weight, factor count).  Finite TopPair checks pass on the
          eight-label roots, the eleven-factor root, F1 multiplicities 2..6
          and the runs of length 5..12.
        - Open: the same identity with individually weighted channels
          `beta_(uvc eta) >= 0`.
      - FM-STR10 (luna_max_neptune; `fm39/str10_support_graph_hall.py`,
        rerun: PASS, about 4 min).  Channel support graph of the FM-STR8b
        correction `delta`, and a Hall test.
        - On 16 exact graphs (the 10 corrected screen profiles and six more
          with labels `<= 5`), `rank delta` equals the multiplicity-weighted
          matching number of the support graph `G` (odd channels with
          capacity `-f_A f_B`, even channels with `f_A f_B`), and both equal
          `dim H_odd(d_0)`.  So Hall holds on all of them.
        - `G` is broad but not complete.  It is not local: `(+1)(-2)^2(+3)^5`
          has an edge `(1,6) -> (7,0)` at distance 6, so no rule of radius
          `<= 5` works.  The unsigned labels do not determine it:
          `(+1,-2,-2,+3,-4,-4)` and `(+1,+2,+2,+3,-4,-4)` have the same
          source capacities and different neighbourhoods.
        - Without an edge rule, the large lists get only the total
          dimension check, which is FM3 itself.  Runs: `O/E` is 6% at k = 6,
          8% at k = 10, 14% at k = 11, 21% at k = 15 and 31% at k = 19
          (`E = 5.43e15`).
          `H_odd(d_0) = 0` for the runs with `k <= 20` not listed.  Random
          lists with `W <= 80`: 25 of 100 have `O > 0`, smallest margin 58.
        - Open (knobs: arbitrary factor count and label size): Hall on the
          actual support graphs, and `rank delta` = matching number.
      - FM-STR8d (astra_max_vulcan; `fm39/str8d_complementary_pair_obstruction.py`,
        rerun: PASS completed STR8d verifier, about 1 min).  Coordinate
        triangular minors are impossible in general.  This does not refute
        the correction itself.
        - Proposition 1 (knob: one even vector per odd vector, scalar
          triangularity): two retained odd vectors with complementary cuts
          `S`, `S^c` are the same vector `v` after forgetting colours.  Their
          columns are `kappa(S cap T) <t,v>` and `kappa(S^c cap T) <t,v>`, so
          they have identical support, and no triangular minor contains
          both.  Example `(-1,-1,-1,-2,+3,+4)`, `(E,O) = (22,2)`.  A dense 2x2
          minor, determinant 2150, still gives `H(d) = (20,0)`.
        - Proposition 2: the obstruction holds for every `(-1)^3(-2)^3(-a)^3
          (-(a+1))^3`, `a >= 4` (for `a = 4`, `(E,O) = (60,756, 1,308)`).
        - Box audit of the 531 corrected profiles: 263 have the
          complementary-pair obstruction; 165 have a triangular minor
          (modulo 1,000,003); singleton peeling stalls on 103.  Allowing 2x2
          blocks covers 71 more.  `(-1)^3(-2)^3(-3)^3(-4)^3`, `(E,O) =
          (16,650, 224)`, still stalls.  All runs k = 5..9 have the
          obstruction (k = 9: 2,988 of 7,664 odd columns lie in
          complementary pairs).
        - Open: `rank(P_even J P_odd) = dim H_odd(d_0)` on the family of
          Proposition 2, and in general.
        - Main-agent note: FM3 needs only one square-zero correction.  For
          any odd-to-even `J`, `delta = PJP` is square-zero and commutes
          with `d_0`.  So generic position weights suffice.  FM-STR8e
          follows this up.
      - FM-SEC178 (luna_max_venus; stopped by the main agent at k = 43 as
        grinding).  The (CB) relative margin on the runs `(-1,...,-k) + (sigma
        p)` rises from 0.98316 (k = 17) to 0.99127 (k = 43), exact rationals.
        The TopPair child/parent ratio at k = 41, p = 43 is about `7.2e-4`.  The
        largest flip there is about `-1.44e-4` times `g_p(B)`.
      - FM-STR11 (luna_max_venus; `fm39/str11_compatible_split_census.py`,
        rerun: ALL ASSERTIONS PASS, 8 min; main-agent check
        `fm39/str11_noflip_split_kill.py`, rerun: FM-STR11 NO-FLIP SPLIT CHECK
        PASS, 10 min).  Channel-compatible splits, `Phi = sum_alpha f_A f_B`.
        - Census: all 19,018 pair-free profiles with labels `<= 5` and length
          2..10 have a compatible nontrivial split; so do 300 random
          profiles.
        - Main agent: this is FM3 itself.  A one-factor peel `(eps n)` is
          compatible exactly when `f_B(n,0) >= 0`, and `Phi(Lambda) = 2
          f_B(n,0)` (the `g_p` form).  A two-factor peel `{u,v}` is
          compatible exactly when `D_uv >= 0` and the fused children are
          nonnegative, i.e. a flip descent.
        - Kill (knob: channel compatibility of a split with both halves of
          size `>= 2`): on ALL 5,430 no-flip lists of the `W <= 40` census, no
          such split is compatible.  So compatible splits give nothing
          beyond flip descent on the hard lists.  The route is closed.
      - FM-STR7e (astra_max_minerva; `fm39/str7e_interior_budget.py`, rerun:
        FM-STR7e PASS, about 1 min; reads `sec166_census_w40_noflip.log.gz`).
        NEW LEAD: an interior height budget gives an explicit injection on
        every tested no-flip list, multiplicity-free runs included.
        - Lemma 1 (knob: grading by the final channel height only): the
          endpoint grading is equivalent to `Phi >= 0` and adds nothing.
        - Lemma 2: for two minus factors, the `Sp(4)` prefixes are
          simultaneous label lowerings, `2<J_(a+b-2-2j) C, K_a K_b> =
          Phi(C, -(a-j), -(b-j))`.  Individual layers can be negative on a
          no-flip list: `(-1,...,-7,-12)`, pair `(-7,-12)`, top layer -4.
        - Proposition 3 (the budget).  For a pair `(u,v)`, `C = Lambda - u -
          v`, cut `C` by the deterministic interior rule: sort by decreasing
          label, assign each factor to the lighter block (ties to A),
          then swap so `weight(A) <= weight(B)`.  With `G_pure = U_a(x)U_b(x) +
          eps_u eps_v U_a(y)U_b(y)` and `G_mix = eps_v U_a(x)U_b(y) + eps_u
          U_a(y)U_b(x)`, put `P_t = sum_(r+s=t) r_A(r,s)[F_B G_pure]_(r,s)`
          and `M_t` likewise with `G_mix`.  If `B_T = sum_(t <= T)[P_t +
          min(M_t, 0)] >= 0` for every `T`, a downward height matching gives a
          sign-reversing involution with `Phi(Lambda)` positive fixed points.
          It may change channels and move many factors.  Since `B_infinity
          <= Phi`, the budget is stronger than FM3 and needs no induction.
        - Proposition 4: no failure in 3,541,270 prefix checks, for EVERY
          pair, on the 5,430 no-flip lists at `W <= 40` (148,972 pairs).
          Also on the runs `k <= 20` with `k+1 <= p <= k+10` (68 lists),
          all 32 signings of F1, and the sign minimum of 1..8 (explicit
          matching: 936 pure plus 20 mixed fixed points, total 956).
        - Theorem 5 (proved): if a cut has `F_A = d^r G`, `F_B = d^(2-r) H`
          with `G, H` genuine `Sp(4)` characters, `r` in {0,1,2}, every
          interior prefix is nonnegative.  This gives a uniform class of
          multiplicity-free words with unboundedly many factors.  It does
          not reach the runs: `d^(2r-2) prod K_n` is not genuine for `r >= 2`,
          and `Phi(d^4 chi_(1,1)) = -6`.
        - Open: for every no-flip residual, some pair (TopPair passed every
          test) has `B_T >= 0` for every `T`.
      - FM-CHK112 (luna_max_eris, fresh code; `fm39/chk112_str9_independent_check.py`,
        rerun: PASS, 2.5 min).  ACCEPT FM-STR9, items 1-4.
        - Proposition 1: matches direct signed tables in 796 random cases.
        - Theorem 2 proof read line by line (NS survival, orthogonality for
          `n != m`, overlap `1/(q+1)`, the 2x2 Gram bound, Bertrand).  Energy
          bound tested exactly on 6,158 TC+NS lists with labels `<= 12` and
          length `<= 8`; 3,121 have nonzero `H_odd(d_0)`, and all pass.
        - Corollary 3 checked for `(3,5,7), M = 18` and `(3,5,9), M = 20`, l =
          2..6.  Proposition 4 dimensions reproduced.
        - So FM3 holds, through the complex, on the TC+NS two-minus class,
          which contains lists with unboundedly many factors.
      - FM-STR8e (vulcan thread, now on Luna; `fm39/str8e_generic_weights.py`,
        rerun: PASS all checks).  Generic weights for the FM-STR8b
        correction.
        - Accepted reduction: for any odd-to-even `J`, `delta = PJP` is a
          differential, so FM3 needs only that SOME weight vector gives
          `rank delta = O`, i.e. generic rank.  The weights `w_i =
          t^((O+1)^i)` keep every nonzero minor polynomial nonzero.
        - Rank `O` at random weights and at `w_i = i+2` on all 268 box
          jobs (the 103 stalled profiles included), the sign minimum (216),
          `(-1)^3(-2)^3(-3)^3(-4)^3` (224), `Lambda_4` (1,308) and the runs
          k = 5..9.  F1 was out of memory range.
        - The separated-scale leading matrix (`e_i = 2^i`) has full rank on
          258 of 268 box jobs.  It fails, for example, on
          `(-1^5,-2,-3^3,-4)` (rank 448 of 490), because tied maximum-weight
          assignments cancel.  32 exponent orders reach 222 of 268.
        - Two-column minor for `Lambda_a`: `H(a)^2 (prod_S w - prod_(S^c) w)`,
          `H(a) > 0`.  Cube condition: when the retained overlap matrix is
          `c prod_(S cap T) w_i`, the leading matrix is a permutation.
        - Knob (support-only matching): a retained 2x2 block with
          proportional rows `(5 kappa(S), 5 kappa(S^c))`, `(5/2 kappa(S),
          5/2 kappa(S^c))` has full support and determinant zero for every
          `w`.  So rank = matching number cannot follow from support alone.
      - FM-STR9c (ceres thread, now on Luna; `fm39/str9c_filtration_obstructions.py`,
        rerun: PASS after a one-character repair in the printed verifier, `range(abs(b-n), a+n+1, 2)` ->
        `b+n+1` in the y-fusion loop; the printed version fails its first
        assertion).  Filtrations for an inductive proof of concentration.
        - Obstructions to an even-concentrated `E_1` page:
          - the `d_0`-first page: all 531 corrected box profiles have odd
            `H(d_0)`;
          - cut-size grading: negative graded Euler characteristic on 463
            of 531, e.g. `(-1,-1,-1,-2,+3,+4)`, grade 3 Euler -6;
          - total height: 44 profiles, e.g. `(+1)^5(-2)(+3)(-4)`;
          - Casimir grading: 511 profiles;
          - pair-colour cut: arrows go both ways.
          Later pages can still cancel: the height example has `E_1 =
          (156,6)` and `E_2 = (150,0)`.
        - A separate mixed-colour piece would need nonnegative Euler
          characteristic, but on `(-1,...,-8)` every pair has mixed Euler
          between -130 and -16.
        - Scalar-insertion obstruction (knob: only parent scalar
          inequalities, associativity and nonnegative individual child
          values): a genuine G-character `R` (not shown to come from a
          signed word) has `F_R >= 0` and nonnegative one-fusion children,
          but `F_(R S_2)(1,3) = -2`.  An induction must keep the
          background's factorization and recoupling maps.
      - FM-STR7e main-agent checks (`fm39/str7e_budget_runs_main.py K1 K2`,
        `fm39/str7e_budget_census_stats_main.py`; exact Python integers).
        - Runs `(-1,...,-k) + (sigma p)`, TopPair `(-(k-2), -k)`, `k+1 <= p <=
          k+3` with even weight, k = 5..34: the budget holds in every case.
          `Phi` agrees with FM-STR8c (184, 980, 6,434 for k = 6, 7, 8).
          Profile (k = 16, p = 18): `P_t > 0` at every height.  Where `M_t <
          0` it is at most about 1% of `P_t`; positive `M_t` reach about 12%
          of `P_t` at the top heights.  The minimum of `B_T / Phi` is small
          only because the first heights carry little mass.
        - Census (TopPair, 5,430 no-flip lists, `W <= 40`): the plain prefix
          and the budget hold on all.  The heightwise version fails on 102,
          and `P_t < 0` at some height on 36.  The largest ratio of negative
          mixed mass to positive pure mass is 0.32, at `B =
          (-1,-3,-3,-4,-5,-5,-6)`, `p = -7`.  So the census lists, not the
          runs, are the tight cases for this mechanism.
        - Consumer note: `Phi = B_infinity + sum_t max(M_t, 0) >= B_infinity`.
          For FM3 alone the `T = infinity` case of the plain prefix suffices,
          and that case is FM3.  The prefix conditions are what the
          involution and the `Sp(4)` truncation method need.
      - FM-STR10b (luna_max_neptune; `fm39/str10b_edge_rule.py`, rerun: PASS,
        about 6 min).  Edge rule for the support graph.  The line pauses.
        - The exact edge rule is the retained-copy overlap formula itself:
          `alpha -> beta` iff some overlap `sum_e c_(T,j)(e) c_(S,i)(e) / prod
          binom(n_l, e_l)` is nonzero.  The four-block CG filter
          `F(S cap T), F((S u T)^c), F(S - T), F(T - S)` sharing a spin is
          necessary, not sufficient: in `(-1,-1,-1,-2,-2,-2,+3)`, 110 of 262
          filtered pairs have zero overlap.
        - On 46 small profiles, the `k = 6` run (44 edges, radius 10) and
          `(+1)^6(-2)(+4)(-6)`: rank = capacitated max-flow = `dim H_odd(d_0)`,
          and one-vertex Hall holds.
        - Odd fraction on `(+1)^(2t)(-2)(+4)(-6)`: 17.8% at t = 12.
        - Knob: computational growth.  Hall at scale needs the full overlap
          matrix, so it is no cheaper than the rank.  No counterexample.
      - FM-STR9d (ceres thread on Luna; `fm39/str9d_two_minus_prefix.py`,
        rerun: PASS FM-STR9d exact verifier, 49 s).  The interior prefix for
        two minus factors with a genuine plus background.
        - Take the two minus factors as the pair and cut the plus factors
          by the interior rule.  Then `F_A = G_A` and `F_B' = d^2 G_B K_n
          K_m`.  If `G_A`, `G_B` (products of `S_a`) are genuine `Sp(4)`
          characters, FM-STR7e Theorem 5 (r = 0) gives every plain prefix
          `>= 0`.  This holds for `(-n)(-m)(+1)^a`, all `n, m, a`, since
          `S_1 = h_1`.  The scalar FM3 values of that family already follow
          from the `Sp(4)` form `2<K_n K_m, h_1^a>`.  The new content is the
          prefix mechanism.
        - Obstruction to the blockwise condition (knob: each block
          genuine): `S_2 = chi_(2,0) - chi_(1,1) + chi_(0,0)` is virtual.
          No negative prefix in 1,134 fundamental-background cases or in
          16,170 bounded arbitrary-plus cases.
        - Open: the plain prefix for every two-minus list, without
          genuine blocks.
      - FM-STR7f (minerva thread on Luna; `fm39/str7f_run_budget_and_block_obstruction.py`,
        rerun: FM-STR7f verifier PASS, 16 s).  The runs budget: no uniform
        proof yet, plus an obstruction to the genuine-block method.
        - Since `P_t + min(M_t,0) = min(P_t, P_t + M_t)`, the stronger LAYERWISE
          condition `P_t >= 0` and `P_t + M_t >= 0` at every height suffices.
          `P_t` is the average of the layers for the pair signs as given and
          flipped.  It holds on the 10 residual runs with k = 7..10 (198
          height checks).  The main-agent profile at k = 16 agrees.
        - Obstruction (knob: at most two extracted `d` factors, and each block
          separately genuine): for `(-1,...,-8)` (residual, `Phi = 980`, all 28
          flips negative), after removing any pair, every block of 3 to 6
          minus labels fails every parity-compatible `d^e G`, `e <= 2`, with
          `G` genuine (210 blocks, 308 quotient checks).  Example: `D_1 D_2 D_3
          / d = chi_(5,0) - chi_(4,1) - chi_(3,2) + 2 chi_(3,0) + 2 chi_(1,0)`.
          The budget itself is positive there.
        - Main agent: `K_n S_a` is genuine exactly when `a <= n` (checked for
          `n <= 8`, `a <= 10`).  Among `S_a S_b` with `a, b <= 6`, only `(1,1)`
          and `(1,2)` are genuine.
        - Open: the runs budget uniformly in `k, p`, or the layerwise
          condition.  A supplier must handle higher powers of `d`, or allow
          cancellation between the two blocks.
      - FM-STR12 (main agent; `fm39/str12_height_prefix_all_splits_main.py`,
        run: FM-STR12 HEIGHT PREFIX ALL SPLITS PASS).  CONJECTURE (HPP):
        height-prefix positivity for EVERY split.  For a signed list with
        an even number of minus factors, every split `A | B` and every `T`,
          `Pi_T(A,B) = sum_(r+s <= T) f_A(r,s) f_B(r,s) >= 0`.
        - `Pi_infinity = Phi(A u B)` (FM3).  `Pi_0 = Phi(A) Phi(B)` when both
          halves are even, else 0.  For `|B| <= 2` HPP reduces to FM3 of the
          fused children (FM-STR7e Lemma 1).  So the content starts at `|A|,
          |B| >= 3`.  In polynomial terms, `Pi_T` pairs the degree-`<= T`
          orthogonal truncations of `F_A` and `F_B` under the product
          semicircle measure.
        - Evidence, no failure:
          - every split of every pair-free list with labels `<= 4`, length
            `<= 8` (79,842 splits);
          - every split with labels `<= 6`, length `<= 7` (227,503);
          - lists with `+-n` pairs allowed, labels `<= 3`, length `<= 8`
            (60,553);
          - every distinct split of 150 sampled no-flip census lists
            (18,090);
          - 800 random splits of the runs `k = 6..18`.
          The interior-cut budget of FM-STR7e is the special case of the
          interior cut.
        - Shape: the weighted height `r + 2s <= T` also passes everywhere.
          Box truncation `r <= R, s <= S`, the `max(r,s)` cutoff and the
          diagonal band `|r - s| <= T` fail (108, 108 and 150 of 22,928
          splits).  A pointwise twisted pairing `sum f_A f_B chi_alpha(k) /
          dim alpha` can be negative.  So HPP is tied to the degree filtration.
        - Consequence: HPP implies every `q`-weighted pairing `sum q^(r+s)
          f_A f_B >= 0`.  That pairing is `<F_A, P_q F_B>` for the Chebyshev
          (free Mehler) Poisson kernel, whose kernel is `(1-q^2)^2 /
          det(1 - q g_1 (x) h_1) det(1 - q g_2 (x) h_2)`.
        - Use: a proof needs only ONE split per list.  HPP says any split
          works, so the prover may choose the most convenient one.
      - ADV-2 (astra_max_vulcan, advisor; `fm39/adv2_constant_weighted_repro.py`,
        rerun exactly).
        - Lemma 1 (decay beyond each label's edge scale):
          `|R_(n,+-)| <= min(1, 4/X)` with `X = (n+1) sqrt(1-z)`; combined
          with FM-MECH49, `|R_(n,+-)| <= (1 + (n+1)^2 (1-z)/4096)^(-1/2)`.
          Main-agent check: 120,000 random high-precision points, no
          violation.
        - Proposition 2 (constant weighted cutoff; independently checked,
          FM-CHK61).  With `T = sum_i (n_i+1)^2/(k+1)^2`, `T >= 2^21` implies
          FM3, with any number of 1's and 2's.  So the plain count needed
          is `H = O(k^2)` instead of `O(k^2 log k)`.
        - Proposition 3: `Theta(k^2)` is optimal for any criterion that
          needs every ray integral to be `>= 0`.  For
          `F(h,n) = E[d^2 hat S_3^(2h) hat S_n]` with `h = m^2`, `n = 10m`,
          the normalized `q = 1` ray tends to
          `1F1(2; 3/2; -5) = -0.02067...`.  Exact witness: `F(4,20) = 406`
          while the ray value is `-3277672/5980345`.  So a constant
          ordinary-count threshold must use the angular average.
        - Lemma 4 / Prop. 5: the Gaussian sphere model is positive for
          every word (sum of squares after `A = (X+Y)/2`, `B = (X-Y)/2`),
          with explicit finite-character errors.  An explicit
          two-dimensional estimate restores positivity along the
          negative-ray family.
        - Distance-3 witnesses
          `F(h, 6h-4) = (2h-1)(4h^2 + 26h + 6)/3` have unbounded `H` and
          values polynomial against exponential scales.  They give no
          lower bound on the FM3 threshold.
        - FM-CHK61 (luna_max_pluto, fresh code): ACCEPT all four items.
          - Lemma 1: 6,480 exact rational profile checks and 400
            angle-formula checks at 150 digits.
          - The product bound, from the chord inequality for `log(1+x)`.
          - Proposition 2: the omitted-factor choice spelled out, and the
            near/far constants recovered exactly.
          - Proposition 3: the exact witness, with the all-`m` error bound
            (7) re-derived.
        - Biggest unknown (all three advisors): a mechanism keeping the
          angular cancellation across mixed label scales, i.e. whether the
          angular average admits a constant ordinary-count `H*`.
      - ADV-1b (astra_max_juno, advisor follow-up; `fm39/adv1b_test1.py`,
        `fm39/adv1b_test2.py`, both rerun exactly): the `H = 1` layer at
        `d >= 6`.
        - Mixed squares survive.  With `c_j^M(xi)` the Krawtchouk-type
          recurrence polynomials (`M = T - 2b`, `xi = a - e`), every requested
          row has an exact rational identity
          `F_d = sum alpha_j c_j^2 + sum beta_(j,nu) (c_j + theta c_(j+2))^2`
          with nonnegative weights.  That covers 26 rows with `T in {13,14}`,
          `b = 3..9`, `d = 6, 7`, plus six rows up to `d = 15`.  The Gram
          matrices are tridiagonal within the even and odd blocks, and the
          identities hold for every real `xi`.  The shift `M = T - b` alone
          fails (negative constant weight at `(36,14,14)`).
        - Exact finite SU(2) midpoint formula (Lemma 2).  If `U` is Haar and
          `V` has density `chi_1(V)^2` with respect to Haar, then `UV` and
          `UV^(-1)` are independent Haar elements.  (Main-agent check: the
          tilted eigenangle density `~ sin^2 theta cos^2 theta` makes `V^2`
          Haar.)  The EVEN form becomes `2^L (||M_0||^2 + <M_0, M_2>)`: a
          square plus one explicit spin-2 coupling.  This is a group-level
          refinement of the torus rotation M0 recorded earlier (whose form
          `f_0^2 + f_0 f_2 - 2 f_1^2` is Lorentzian).
          - For `H = 1` it gives an explicit channel sum over `(l, j)` with
            nonnegative angular weights.  Individual channels can be
            negative, but grouped sums over `l` at fixed `j` (or over `j`
            at fixed `l`) were nonnegative in all 105 evaluated profiles
            through `d = 12`; the four requested rows agree exactly.
        - The pair induction on `e` closes algebraically as a
          two-component transfer system in `(m_0, m_1)` with a `b`-shift.
          Pair positivity alone is equivalent to the target, so an
          invariant is still needed.  `(x-y)^3 = Psi_30 - 3 Psi_21 + 5 Psi_10`
          shows why the antisymmetric cone does not extend directly.
        - KILL (scoped): output reindexing within a common imbalance class;
          separator `-346112/1155` at `(T,b,d) = (13,3,6)`.
        - Recommendations: (1) a recursive (Schur-pivot) construction of the
          tridiagonal mixed Gram matrices; (2) grouped Jacobi positivity
          from the midpoint formula; (3) a transfer-invariant cone for the
          `e`-induction.
      - ADV-1c (astra_max_juno, advisor; `fm39/adv1c_part1.py`
        (`--part small|kernels|large --block k`), `fm39/adv1c_part2.py`, rerun):
        the midpoint identity for the whole cone.
        - Identity (every list): with `C_n = R_n + R_n^*`,
          `S_n = i(R_n - R_n^*)` (plus and minus features), the feature
          tensor `T = tensor of F_i`, `M_0 = int T` and
          `M_2 = int chi_2 T` (Haar), one has
          `EVEN = 2^L (||M_0||^2 + <M_0, M_2>)` (Hilbert-Schmidt).
          Checked exactly on all 6,470 even-minus classes with labels `<= 4`
          and length `<= 8`, on 128 large-label cases, and by 243 quaternion
          kernel integrals.
        - It is a reformulation: with `A = ||M_0||^2` and `B = <M_0, M_2>`,
          FM3 is exactly `B >= -A`.
        - KILLs of the natural strengthenings (exact families):
          - `|B| <= K A` fails for every `K`: `hat S_n hat S_1^n` has
            `M_2 = n M_0`.
          - Norm contraction fails (`hat S_2^3`).
          - Positivity on the feature span fails, already inside the span of
            `hat S_2^3`.
          - Blockwise positivity in total spin fails: `hat S_2^3` has blocks
            `1/8, -1/16, 3/16, 0`.
          - Bare repeated tilting gives no decay.
          The constant 1 is sharp: `D_1^(2m)` has
          `M_2 = -m/(m+2) M_0`.
        - Projector form: `tr[P_0 Pi_eps (P_0 + P_2/3)] >= 0`, with `Pi_eps`
          built from partial time reversals.  A proof must combine
          different total-spin sectors.
        - Biggest unknown (shared with ADV-2, ADV-3): a uniform grouping
          that retains cross-sector cancellation.
      - FM-SEC136 (luna_max_jupiter; `fm39/sec136_gram_part1..3.py`, rerun
        exactly): the tridiagonal mixed Gram family for `H = 1`.
        - With `c_n^M = sum_s rho_(n,s) xi^(n-2s)`, every tridiagonal (even/odd
          block) Gram representation of `F_d(T, u, b)` lies in an affine
          family with `d - 1` free off-diagonal entries `beta_j`.  The diagonal
          `A_m(beta)` follows by a descending recurrence, using
          `E_(m,m) = 1/(m!)^2`.
        - Grid `4 <= T <= 60`, `3 <= b <= 20`, `6 <= d <= 20`, `p >= 7`: the
          diagonal member `beta = 0` has a negative entry in 7,660 of 12,951
          rows (first `(4,8,6)`, `M = -12`).  That row still has an exact mixed
          certificate with positive Schur pivots in both blocks.
        - With the finite slope menu `theta = +-2^s`, 574 of 595 rows are
          covered (`T <= 20`, `b <= 12`, `d <= 10`).  The first uncovered row,
          `(4,12,10)`, has a polynomial with all positive coefficients in
          `u`.  The menu separator there is not an all-slope obstruction.
        - Open: a uniform choice `beta(T, b, d)` with nonnegative Schur
          pivots.
        - FM-SEC137 (luna_max_jupiter; `fm39/sec137_gram_census_repro.py`,
          rerun exactly in 16 min).
          - Bidiagonal LDL in each parity block: `P_r = A_r`,
            `l_j = beta_j/P_j`, `P_(j+2) = A_(j+2) - beta_j^2/P_j`.
          - On the 1,655-row grid (`T <= 30`, `b <= 12`, `6 <= d <= 14`):
            - 1,146 rows are certified exactly (one-start search plus
              rational LDL checks);
            - 509 rows are unresolved by the search, not shown infeasible;
              the first is `(5,12,11)`;
            - no exact obstruction to a tridiagonal PSD member was found.
          - No canonical closed-form `beta` was found.
          - Main-agent check: 543 of the 1,655 rows have all coefficients
            of `F_d` in `u` nonnegative, mostly for `M = T - 2b <= 6`.  Only
            17 of the 509 unresolved rows are among them, leaving 492
            spread over all `M`.
          - Per-row certificates do not finish the layer, so this route is
            parked as evidence.
      - FM-MECH62 (astra_max_juno; `fm39/mech62_grouped_midpoint_kill_repro.py`,
        rerun exactly): KILL of grouped positivity in the midpoint formula.
        - Canonical Jacobi form: with `mu` the semicircle on `[-1,1]`,
          Gegenbauer bases `C_n^(l+1)` and `J_l` = multiplication by `t`,
          each fixed-`(j,l)` group is one selected coordinate of the positive
          operator `J_l^2`, with Legendre weights `w_l >= 0` summing to 1.
          The selected-coordinate form has a 2x2 principal minor with
          negative determinant.
        - Exact counterexamples on residual profiles below the cutoff:
          - the balanced row `(e,a,b,p,d) = (2,2,15,22,6)` has negative
            fixed-`j` groups `j = 9..12` and negative fixed-`l` groups
            `l = 6, 8, ..., 16`, while `F = 504656`;
          - `(23,2,3,19,6)` has fixed-`j` pieces
            `129908066/2925, -55104044146/2925, -235766538154/2925,
            1239529442878/975`, summing to 1171913728.
        - What survives: summing over output labels gives a uniform square
          identity, i.e. a positive average over core labels.  That is
          weaker than the consumer.  Knob: freezing `j` or `l` at a fixed
          output label.
        - So the `H = 1`, `d >= 6` routes killed so far are: raywise
          positivity, a fixed Gaussian margin, diagonal squares, the square
          template, positive reductions to `b = 0`, Newton in `b`, and
          grouped midpoint sums.  Remaining routes: per-row mixed Gram
          (parked) and a transfer invariant for the `e`-induction
          (FM-MECH65).
      - FM-MECH65 (astra_max_juno; `fm39/mech65_transfer_cone_repro.py`,
        rerun exactly): the transfer-cone route for the `H = 1` induction on
        `e`.  Not closed.
        - The exact two-component transfer `(m_0, m_1)` under `x - y`, with
          the `b`-shift, was verified on 275 full U-coefficient vectors.
        - KILLs:
          - fixed-output Hankel positivity: uniform failure on required
            seeds, with a negative leading coefficient at large labels for
            every fixed shift of `b`;
          - coupled block-Hankel: negative quadratic `-2680` on the actual
            orbit into `(6,8,8)`;
          - total positivity of the coefficient kernel;
          - the one-step two-component order cone.
          Knobs: the stated Hankel blocks and order conditions.
        - Positive: the signed-character cone `C_-` is invariant under
          `x - y` and `Z` and gives U-positivity.  But it contains only the
          `a = 0` seeds (already known); `(x+y)` already has coefficient `-2`
          at `(1,1)`.
        - Status of `H = 1`, `d >= 6`: grouped midpoint false; transfer cone
          only for `a = 0`; per-row mixed Gram parked.  What remains
          unexcluded is a construction coupling output labels and `b`.
      - ADV-3 (astra_max_minerva, advisor: gap audit and red team;
        `fm39/adv3_onelabel_allb_search.py`, `fm39/adv3_cp_certificate_1422MM2.py`,
        both rerun exactly).
        - No positivity theorem is affected.  FM-MECH45, 46, 49 and 50
          verifiers were rerun; the FM-MECH47/48 analytic cutoffs replayed.
        - REPAIR (bookkeeping) of the FM-MECH49 Prop. 7 residual filter.
          - (i) LS, LR4 and LLm-S must be applied to the WHOLE
            non-fundamental word, `hat S_2` factors included (LS: each
            `hat S_2` adds 8/3 to `Lambda`; LR4: adds 8 to `L`).  Testing
            only the labels `>= 3` is invalid.  Witness:
            `phi_2(hat S_8 hat S_2^5) = 560` with the nonzero cross term
            `E[d^4 U_2(y) U_8(x) U_2(x)^4] = 33` (both recomputed by the
            main agent).
          - (ii) The displayed necessary conditions are not an exact
            membership test: `phi_5(h_4^4 h_1^2) = 1600` satisfies them, but
            its integrand is a square.
          - (iii) Re-encoding `-2` as `-1, +1` can hide an exclusion: the
            list `(-2,-2,+2,+7,+7)` is covered by FM-MECH39 (distance 3, no
            label 1).  A coverage checker must keep the original list.
        - Targeted exact counterexample searches found no negative and no
          non-support zero:
          - one large label at distances 3..6: 40,431 sign cases;
          - mixed scales with 2..5 core factors: 116,710;
          - labels up to 1,594: 12,090;
          - one-label profiles over `b`: 86,838 U-coefficients, among them
            58,654 residual candidates (`b >= 3`, `n >= 7`, distance `>= 3`),
            least coefficient `[U_7] g_(3,4,3) = 32`;
          - (E) rows: 1,075,304 inequalities.
        - Proposition 1: every list `(1^4, 2^2, M, M+2)`, `M >= 4`, has an
          explicit seven-atom H_AC certificate (nonradial, list-dependent
          autocorrelations), and every FM3 value there is `>= 24`.  B fails
          for `(1^7, 5)` while `f = (1/2) p*p + (5/2) delta_0` (H_AC).
        - Controls.  The common-radius normalization has
          `K^2 (1 - R^2) -> 17/10`, so it cannot give a label-independent
          per-factor loss.  The generating polynomial is not stable
          (`1 + z_1 z_2` vanishes at `z_1 = z_2 = i`).
        - Ranked recommendations: (1) uniform-in-`b` positivity on the actual
          one-label profiles; (2) list-dependent H_AC at q = 0 with
          nonradial factors, starting from Prop. 1 at distance 3;
          (3) exhaustiveness of the (E) criteria for `q >= 3`; (4) a
          label-independent `H*` (low-medium confidence).
      - ADV-1 (astra_max_juno, advisor; `fm39/adv1_distance3_repro.py`, rerun
        exactly).
        - **Lemma 1 (PROVED; main-agent algebra check): the `H = 1` layer at
          list distance 3, for every label and every `b >= 2`.**  With
          `F = [U_n] E_y[(x-y)^e (x+y)^a Z^b]`, `t = a + e`, `u = (a-e)^2` and
          `n = t + 2b - 6 >= 7`:
          `144 F = u^3 + A u^2 + B u + C`, with
          - `A = 30b - 3t - 20`,
          - `B = 9t^2 - 36bt + 180b^2 - 6t - 300b + 64`,
          - `C = 3(3t^3 + 6bt^2 - 18t^2 - 12b^2 t + 48bt + 40b^3 - 120b^2 + 32b)`.
          Positivity:
          - `B = 9(t - 2b - 1/3)^2 + 144(b-2)^2 + 264(b-2) + 15`;
          - `C = 9(t-2)(t^2+16)` at `b = 2`, and `C/3` has a positive form
            at `b = 3 + v`;
          - if `A < 0`, integrality forces `t >= 10b - 6`, and then
            `(4B - A^2)/9 > 0`, so `144F = u[(u + A/2)^2 + (4B - A^2)/4] + C`.
          Main-agent checks: every decomposition recomputed with sympy; the
          identity holds on 263 fresh random cases from
          `extra_label_upositivity_screen.py` (both parities of `e`), on top
          of the advisor's 2,830.
        - Lemma 2 (saturation).  With `n = sum mu_i - 2d`, the q = 0 table
          is `g(S) = m(mu_S) K_(d - w(S)/2)(mu_(S^c))`, where
          `K_j(nu) = [z^j](1-z) prod_i (1 + ... + z^(nu_i))`.  Every
          non-distinguished label `> d` can be clipped to `d + 1` without
          changing the table.  So at fixed distance the arbitrary-label
          dependence disappears, and the background enters as finite-degree
          polynomials in its counts.
        - Gaussian/Bessel limit.  With `X = U + V`, `Y = U - V` (`U`, `V` iid
          isotropic Gaussians in `R^3`), every plus factor becomes a product
          of two cosines and every minus factor minus a product of two
          sines.  The limit is therefore an integral of squares for every
          `H` and every sign pattern.  Finite-label errors look like
          `O(1/a)` (diagnostics only).
        - KILLs.
          - Termwise cyclotomic-row positivity: the `b = 2` witness
            `phi_2(hat S_7 h_1^5 hat S_2^2) = 44` has row contributions
            4, 8, -4.
          - Raw total positivity of the kernel.
          - The canonical invariant-space operator `sum_S eps_S P_S`: it is
            `diag(3, -1)` for four fundamental labels with two minus signs.
        - Ranked recommendations: (1) fixed-distance scalar polynomials
          with a recurrence in `d`; (2) a summed transfer invariant on the
          actual Krawtchouk profiles; (3) the Gaussian squares with relative
          finite-label corrections.
      - FM-SEC127 (luna_max_mars; `fm39/sec127_q2plus_audit_repro.py`, rerun
        exactly): the q = 2 plus sign of (E), uniformly in the label.
        - With `x = c_j`, `y = c_(j+1)` and the Krawtchouk recurrence, the
          slack `E+ = D_j - D_(j+3) + W` is a quadratic form
          `A x^2 + B xy + C y^2` whose coefficients depend only on
          `(N, a - e, j)`.  `W >= 0` gives `E+ >= 0` by Theorem OL.
        - Census `a, e <= 40`: 33,200 rows, no negative slack.  All 20,715
          rows with `W < 0` have a positive definite quadratic.
        - Main-agent extension (`fm39/sec127_pd_dichotomy_check.py`,
          `a, e <= 70`, 175,175 rows).  The dichotomy "PD or `W >= 0`" holds
          everywhere: 14,420 rows are not PD, and none of them has `W < 0`.
          `A > 0` in every row.  Non-PD rows need `|a-e|/N >= 0.32` and
          `2j - N >= 6` (label `>= 8`), apart from 6 rows at `2j = N`.  This
          is the lopsided region where Theorem OL's outer-region ratio
          bounds apply.
        - So a uniform proof splits into (i) `A > 0` and an explicit
          description of the region `B^2 >= 4AC`, and (ii) `W >= 0` there,
          from the actual ratio `c_(j+1)/c_j`.  Delegated as FM-SEC129.
        - FM-SEC129 (luna_max_mars; `fm39/sec129_q2plus_regions_repro.py`,
          rerun exactly).  Fold to `a >= e`; put `d = a - e`, `n = N - j`,
          `kappa = 2j - N`, `t = e`.
          - (i) PROVED: `A = (kappa+2) alpha/H > 0` with explicit
            positive `alpha`, `H = (n+1)(j+2)^2(j+3)^2(j+4)`.
          - (ii) The non-PD region is
            `R = {Delta_num = d^2 beta^2 - 4(kappa+2) alpha gamma >= 0}`, with
            `beta`, `gamma` explicit (in the reproducer).
          - `W = x^2 Omega(r)`, `r = c_(j+1)/c_j`, with `Omega` quadratic in
            `r` (coefficients `omega_0`, `omega_1`, `omega_2` explicit).
          - (iii) Box `a, e <= 70`: all 14,420 rows of `R` lie in Theorem
            OL's outer region `d^2 >= 4(n-1)(j+2)` (14,394) or its central
            region `N >= 3t(kappa+1)^2` (26 central-only).  `Omega >= 0` on
            the OL outer ratio interval `0 < r <= R_OL` in every outer row.
          - Open, in closed form: `R_miss = R minus (outer union central)`
            is empty; `Omega >= 0` on `(0, R_OL]` over `R` intersect outer;
            `Omega >= 0` on the central ratio ranges.  Delegated as FM-SEC131.
        - FM-SEC131 (luna_max_mars; `fm39/sec131_n1_slice_repro.py`).
          - The `n = 1` outer slice is proved:
            `Omega = (1 - d r)(omega_0 - omega_2 r/d)`, and
            `Delta_num >= 0` forces `d^2 >= 3k + 16`.
          - The `k = 0`, odd-`t` zero-coefficient boundary is handled.
          - Open: the coverage step (a), and (b), (c) for `n >= 2`.
        - Main-agent scan (`fm39/sec131_band_scan.py`).
          - The band `R` minus the outer region is thin: the smallest `d`
            with `Delta_num >= 0` is `>= 0.9475 d_out`.
          - It occurs only at `k <= 8` for `n < 400` (9,619 points), with
            `t` growing slowly (up to 15 at `n, k <= 300`).
          - Every band point satisfies the central condition, with
            `N/(3t(k+1)^2) >= 5` in the sample.
          - So (a) should split into two steps: `R` lies inside the outer
            region for `k >= K_0`, and the central condition holds on the
            band for `k < K_0`.
        - FM-SEC132 (luna_max_mars; `fm39/sec132_t3_family_repro.py`, rerun).
          - KILL of the main-agent split (a1): for `t = 3`,
            `d = 2n + k - 6`, one has `[n^10] Delta_num = 144`, and the band
            persists at every `k` (e.g. `n = 1283`, `k = 9`).  These rows
            are central, so the coverage claim (a) survives.  Knob: my
            proposed intermediate, not the consumer.
          - `omega_2 = k(j+4)(D - T_2) >= 0` on the whole outer region,
            since `D_out - T_2 = 3(n-1)(n+k+2)`.
          - Outer endpoint reduction: with `omega_2 >= 0`, it suffices to
            show `Q'(R_OL) <= 0` and `Q(R_OL) >= 0`, both in explicit
            `sqrt S` form.
        - FM-SEC133 (luna_max_mars; `fm39/sec133_gap_repro.py`).
          - Parametrize by `t`, with `d = 2n + k - 2t`.  Then
            `D_out - D = 8nt + 4n + 4kt - 4k - k^2 - 4t^2 - 8`.
          - A second central band family: `t = 3`, `n = 16q^2`, `k = q`.
          - Coverage (a) is equivalent to `Delta_num(n,k,2n+k-2t) < 0` on
            the gap
            `(k^2 + 4k + 4t^2 + 8 - 4kt)/(8t+4) < n < (3t(k+1)^2 - k)/2`,
            `t >= 3`.  This is screened for `3 <= t <= 7`, `k <= 15`;
            it is not proved.
          - Paused after three rounds (FM-SEC131..133) without closing a
            uniform step.  It will be re-prioritized once FM-MECH49 states
            the full-cone residual.
      - FM-MECH46 (astra_max_ceres; `fm39/mech46_insertion_repro.py`, rerun
        exactly): the uniform one-extra-label statement, the b-insertion.
        - Exact recurrence (Prop. 1): with `g_b = g_(k,a,b)`,
          `D = k + a + 2b`, `R_D = (D + 4 - x d_x)^(-1)` (acting as
          `x^j -> x^j/(D+4-j)`):
          `g_(b+1) = (x^2+2) g_b + 4 R_D[ b (x^2+2) g_(b-1) - (b+3) g_b ]`
          (290 exact checks).
        - Lemma 2: `R_D` preserves U-positivity (descending induction on
          the U-coefficients).
        - KILL (Theorem 3): no linear insertion `V_L -> V_(L+2)` sending every
          `g_(k,L-k,0)` to `g_(k,L-k,1)` preserves the whole U-positive cone.
          The map is forced to be `T_L = (x^2+2) - 12 R_L`, and
          `T_3 U_3 = U_5 + U_3 - U_1`; for `L >= 4`,
          `[U_(L-4)] T_L U_L = -(L-3)(L+2)/4`, at distance 3 for `L = 12`.
          Knob: positivity on the whole cone, stronger than the consumer
          (the actual children `g_(2,1,1)`, `g_(3,0,1)` are positive).
        - Prop. 4: for two extra labels on the background
          `F = d^k s^a Z^b`, positivity for both compatible signs is exactly
          `R_(n,m) = sum_(j in CG(n,m)) c_(j,0) >= |c_(n,m)|`,
          `c_(n,m) = E[F U_n(x) U_m(y)]` (6,480 exact checks, least slack 1).
          At `b = 0` this is (E) of FM-MECH7.
        - Open: the insertion inequality for the actual profiles,
          `4(b+3) [U_n] R_D g_b <= [U_n](x^2+2) g_b
          + 4b [U_n] R_D((x^2+2) g_(b-1))`, relevant for `n >= 7`,
          `n <= D - 4`.
        - FM-SEC125 (luna_max_mars; `fm39/sec125_label4_moments_repro.py`,
          rerun exactly): the label-4 inequalities of FM-MECH44 in joint
          moments `mu_(m,k)`.  Proved for `b = 0, 1` by explicit
          coefficientwise-positive numerators; 4,851 exact cases with
          `A, E, b <= 20` are all `>= 0`.  Superseded by FM-MECH45 (all `b`);
          kept as an independent corroboration.  The block expansion
          alternates, so positivity only appears after summing.
        - Main-agent screen (`fm39/sector123_moment_screen.py`): the
          {1,2,3} sector.  `E[s^A d^E Z^b (Z+P)^al (Z-P)^ga] >= 0` in all
          7,182 cases with `A, E` even and `A + E + 2b + 2al + 2ga <= 26`.
          This includes 4,139 cases outside the FM3 constraints `al <= E`,
          `ga <= A`, i.e. GKS2*-type words with extra `h_2` factors.
          Delegated as FM-SEC126 (screens, mechanisms) and FM-MECH47 (proof).
          - On a ray, `Z +- P = (1 +- c) f +- 2c`, equal to
            `2((1 +- q)^2 + q^2) > 0` at `f = B`.  So FM-MECH45 Prop. 3/7
            give each fixed `(al, ga)` outside a finite box.  At `q = 0`
            the criterion needs about `N + b + 3 >= 4(al + ga)`, so it
            misses the words with many 3's.  Pure `(3^L)` is B by FM-MECH31.
            What is open is the range in between, uniformly in the number
            of 3's.
        - FM-SEC126 (luna_max_mercury; `fm39/sec126_sector123_screen_repro.py`,
          rerun exactly in 10 s).
          - Exact recurrences in `(b, al, ga)` over the joint moments
            `mu(m,k)`.
          - Census with total degree `<= 60`: 324,632 tuples (145,167
            inside the FM3 constraints, 179,465 outside), no negative value.
          - Eight zeros, among them `(A,E,b,al,ga) = (0,0,0,1,2)`,
            `(0,0,0,2,1)` and `(0,0,1,1,1)`.  The minimum of value/`mu` is
            1/5, at `(4,0,0,0,1)`.
          - KILL: the moment-positive cone
            `{G : E[s^A d^E G] >= 0 for all even A, E}` is not closed under
            multiplication by `Z +- P`.  With `H = Z + P`, `G = (H-2)^2`:
            `E[H G] = 8 - 12 + 0 = -4`.  Knob: the whole cone; the target
            monomials stay positive.
          - Main-agent test: the naive induction hypothesis "`F/mu` is
            nondecreasing in `m` and `k`" fails (2,625 of 5,880 pairs).  So
            a ratio induction needs direction-dependent weights.
        - Main-agent reductions and sharpness
          (`fm39/sector123_quadratic_bstar.py`).
          - `Z = ((Z+P) + (Z-P))/2`, so in the unconstrained form the
            `Z^b` factors are free: the target is
            `E[s^A d^E (Z+P)^al (Z-P)^ga] >= 0`.
          - Multi-affinity in the cross coefficients makes this equivalent
            to `E[s^A d^E prod_i (x^2 + 2 b_i xy + y^2 - 2)] >= 0` for all
            `b_i in [-1/2, 1/2]`.
          - The critical half-width `B*` (largest `B` with every sign
            vertex `>= 0`, `A, E <= 8`) is:
            - exactly `1/2` for 3 factors at `A = E = 0`; this is the zero
              `E[(Z+P)(Z-P)^2] = 0`;
            - `0.6124` for 2 factors, and 0.5774, 0.5968, 0.6219 for
              5, 7, 9 factors;
            - larger for even factor counts.
          - So `b = +-1/2` (the label-3 quadratics) sit exactly on the
            boundary, and any proof of the unconstrained form must be
            sharp at that zero.
      - FM-SEC123 (luna_max_venus; `fm39/sec123_admissible_sets_repro.py`):
        admissible label sets.
        - Among {2,4}, {2,3}, {3}, {4}, {2,4,6}, {2,3,4}, the even labels and
          the odd labels `>= 3`, only {2} is proved (by FM-MECH41).
        - Law-level condition: if `sgn(X)` is a character, the `U_(2k+1)`
          coefficient there is `8(-1)^k (k+1)/(pi(2k+1)(2k+3))`.  This is
          negative for `n = 3 (mod 4)`, so a character sign is impossible
          for any set containing 3 or 7.
        - The quantile encoding has a negative `U_4` coefficient at
          `sigma_1`, `sin(4 theta_0)/pi - 2 sin(6 theta_0)/(3 pi)`.
        - B checks: `(2^a, 4^b)` for `a + b <= 6` and `(4^M)` for
          `M <= 9` all pass (30 exact certificates).
        - Open: `{2,4}`, which needs `b >= 0` and `b*b - b >= 0` off 0 for
          the coefficients of `Y = U_2(X)` (with the disjoint-support
          condition, `b*b = b` on `supp b`); and a non-character-sign
          realization for 3.
      - FM-MECH43 (astra_max_ceres; `fm39/mech43_extra_label_repro.py`,
        rerun exactly): the {1,2} sector plus one extra label.
        - B certificates for eight lists `(1^a, 2^b, n)` with `n = 3, 4`.
        - IBP lemma: with `s = x+y`, `d = x-y`, `Z = x^2+y^2-2`, `P = xy`,
          `M_b = E[s^A d^E Z^b]` and `Q_b = E[s^A d^E P Z^b]`,
          `(A+E+2b+6) Q_b = 4(A-E) M_b + 12 b Q_(b-1)`.  This follows from
          `E[(4-x^2) d_x F] = 3 E[x F]` applied to
          `V = y(4-x^2) d_x + x(4-y^2) d_y`, so `Q_b` is an explicit
          positive combination of the `M_j`.
        - Since `h_2 = Z + P` (FM-MECH44 correction: no `+1`) and
          `hat S_3 = s(Z - P)`, the whole
          extra-label-3 sector `phi_r(h_2 h_1^a hat S_2^b)`,
          `phi_r(hat S_3 h_1^a hat S_2^b)` reduces to the single inequality
          `M_(b+1) >= Q_b` for even `A > E >= 2`.  It passes 637 exact
          checks; the boundary `E = 0` and the case `A <= E` are
          immediate.
        - KILL: deriving it from positivity of the `M_j` alone.  With
          `M_j = 2^j` (the moments of `delta_2`) it fails by `-512/105`, so
          the actual joint law is needed.
        - FM-SEC124 (luna_max_mars; `fm39/sec124_extra_label_sum_repro.py`).
          Expanding in the joint moments gives inner blocks
          `H_(m,k) = 2(m^2 - m + 5k^2 + 7k - 2mk) mu_(m,k)/((m+k+2)(m+k+3))
          >= 0`, but the outer sum over `Z^b` alternates (termwise failure at
          `(4,2,1)`; the true value is `+8`).  945 exact checks
          (`A, E <= 20`, `b <= 20`), all positive.
        - Main-agent restatement: `-P = (d^2 - s^2)/4`, so
          `M_(b+1) - Q_b = E[s^A d^E Z^(b+1)] + (1/4) E[s^A d^(E+2) Z^b]
          - (1/4) E[s^(A+2) d^E Z^b]`.  All three are {1,2}-sector EVEN
          values, each `>= 0` by FM-MECH41.  So the gap is a comparison
          between {1,2} values: two plus-1's cost at most one plus-2 plus
          two minus-1's.
      - FM-SEC121 (luna_max_venus; `fm39/sec121_B_small_labels_repro.py`,
        rerun exactly):
        B census for small labels.
        - Every even-total list with labels `<= 4` at `L = 7, 8, 9` (60, 85
          and 110 lists) has an exact subgroup-orbit B certificate, with up
          to 25 terms.
        - Every list with labels `<= 3` at `L = 9` is certified, and five
          samples at `L = 10`.  No Fourier-negative entry.
        - So B for all lists with labels `<= 4` is consistent with the data.
          By FM-MECH42 it cannot come from one common realization, so it
          would need list-dependent structure.
        - First unscreened case: `(1^6, 3^4)`.
      - **FM-MECH40 (astra_max_ceres;
        `fm39/mech40_insertion_obstruction_repro.py`, all assertions
        pass): no parent-preserving insertion.**
        - Setting: insert a label `a` into `mu` with `n` fixed.  The child
          slices are `A = sum_(j in CG(a,n)) A_j` and `C`, and the parent is
          `P`.
        - A construction `P = sum c Q(u)`, `A = sum c [Q(u) + Q(v)]`,
          `C = 2 sum c u*v` must satisfy, for every character `T` and
          q-degree `k`, `C_k(T)^2 <= 4 P_k(T) (A_k(T) - P_k(T))`.  This is
          a Gram matrix bound, so it holds even for signed factors and
          arbitrary partners.
        - KILL: the bound fails at `(2^8; a = n = 2)`, i.e.
          `(2^9) -> (2^10)`, at `q^4` (gap `-506,106,194,944`), and at
          `q = 0` for `(2^15)`.  Asymptotically, tilting by `(Z+Z')^m`
          (`Z = U_2`) gives ratio `-1`.
        - Both tables are H_AC_q (83 control certificates).  Reorganizing
          the fusion channels (each `(A_j, C_j)` H_AC_q with
          `sum_j C_j = C`) passes the census: 47 insertion instances, 760
          certificates, and a 24-term allocation at the obstructing `q^4`.
        - Knob: preservation of the parent autocorrelation.  An induction
          on labels must re-factorize the old half; no
          'keep the parent, add partners' rule reaches the full cone.
      - FM-SEC115 (luna_max_mercury; `fm39/sec115_hACq_d2_repro.py`, rerun
        exactly): an independent proof of H_AC_q at distance two.
        - Explicit table: `D_mu` at the origin, `alpha` on type 11,
          `1+q` on types 22 and 112, `2+q` on type 1111.
        - Factors `E, A, B, P_j = A + (3v/2) delta_(e_j)`,
          `Q_j = A + 3v delta_(e_j)`, `R_jk`.
        - The residual `rho(q)` is coefficientwise nonnegative by an
          append induction: `J_(nu a) = J_nu + [S][a]`,
          `D_(nu a) = D_nu + J_nu [S-2][a] + [S][S-1][a choose 2]`.
        - Census: 3,080 ordered lists, 496,584 entries.  This
          corroborates FM-MECH38.
      - **FM-MECH39 (astra_max_ceres; `fm39/mech39_hACq_general_repro.py`,
        both parts rerun exactly): H_AC_q uniformly for all labels
        `>= d`, and at `d = 3` with at most one label 1.**
        - Lemma (general contraction formula): the coefficient
          `K_d(mu) = [H_n] prod H_(mu_i)` is `sum over partial matchings M`
          of the strands (`d` edges, none inside a block) of
          `q^(cr(M) + cov(M))`.  Here `cov` counts unmatched strands
          lying inside an edge.  Equivalently it is a product of
          q-binomials and `[k]_q!` along the blocks.  So the free-strand
          reading is exact: a summand of `g_q(S)` is a complete matching
          on `S` and a partial matching on `S^c`, with `2d` contracted
          strands in total.
        - Theorem (uniform in `d`): if every non-distinguished label is
          `>= d`, then
          `g_q = ([d]_q!/2) Q(E) + (K_d - t [d]_q!/2) delta_0`, with `E` the
          labels equal to `d`, and `K_d >= (L-1)[d]_q!`.  So H_AC_q, hence
          FM3, holds on that region at every distance.
        - Theorem (`d = 3`, at most one label 1): explicit factors
          `B, V, T` (triples of 2's), `B + delta_e * V`.  The origin
          residual is nonnegative by the recurrence (3) on `2^v` and the
          append/grow identities (4), (5).
        - KILL: origin-only transfers between averaged trajectories.  At
          `(1^4, 2^3)`, `d = 3`, `q^0` the origin budget is 13 against
          `4321/315`.  Four such lists occur in the `L <= 7` census.
        - `(1^4, 2^2)` has a q-positive repair.
        - Open: mixed `d = 3` with two or more labels 1, and general `d`
          with labels below `d`.
        - FM-CHK52 (luna_max_eris, fresh code).
          - ACCEPT: Lemma 1 (24 random lists, `d <= 4`).
          - ACCEPT: Theorem 2, on 723 lists with `d = 1..6` and 15,294
            entries.
          - ACCEPT: Theorem 3, on 808 lists, including the exception
            `(1,2,2,2)` via `(1,1,2,2)`.
          - REPAIR, Propositions 4 and 5: the four transfer failures and
            the `(1^4,2^2)` repair are verified, but the checker's own
            census did not finish.  The main agent's rerun of the
            reproducer's `--part transfers` did complete through
            `L = 7` (`[20, 347, 4]`), so these four are the only failures
            in that range.
      - **FM-SEC111 (luna_max_pluto; `fm39/sec111_bernstein_hAC_repro.py`,
        rerun exactly): conjecture BH, the Bernstein tables along `q = 0`
        are H_AC.**
        - Write `f_(0,s) = sum_j C(d,j) s^j (1-s)^(d-j) g_j`.
        - Census: all 545 ordered even-total lists of length `<= 6`, labels
          `<= 3`: 3,757 tables `g_j`.  All have exact subgroup-sum
          certificates, with no negative entry and no negative Fourier
          value.
        - `(1^4, 2^3)`: every `g_j` is a subgroup sum.  `(1^8,6)`: `g_0..g_5`
          are subgroup sums, and `g_6` (the FM3 table) is outside B (exact
          separator over all 417,199 subspaces) but H_AC,
          `3 delta_0 + (1/2) p_1*p_1`.
        - The increments `g_(j+1) - g_j` are not H_AC (at `(1^4)`,
          `g_1 - g_0 = delta_5`).
        - The step is a subgroup merge: for `H cap K = {0}`,
          `1_H + 1_K + 1_((H\0)+(K\0)) = 1_(H+K) + delta_0` (the Pluecker
          triple).  It gives `g_0 -> g_1` at `(1^4)` and, with weights, at
          `(1^8,6)`.
        - Open: a uniform merge rule, including the final step where B
          gives way to H_AC.
      - FM-SEC118 (luna_max_venus; `fm39/sec118_BHq_repro.py`; the
        `smallcert` and `boundary9` stages rerun exactly): conjecture BH_q.
        - Statement: the Bernstein tables in `s` of the braided pair table
          `f_(q,s)` are H_AC_q, with one q-independent dictionary.  This
          contains FM3, H_AC, H_AC_q and BH.
        - Census: all 545 ordered even-total lists of length `<= 6`,
          labels `<= 3`, and 90 representatives under rotation and
          reversal: 3,757 Bernstein tables, 63,986 q-coefficient vectors.
          No negative entry or Fourier value, and every nonzero table is a
          nonnegative sum of subgroup indicators.
        - `(1^4, 2^3)`: subgroup sums throughout.  `(1^8,6)`: every table
          is a subgroup sum except the last one at `q^0` and `q^1`, which
          need the sphere: `3 delta_0 + (1/2) p_1*p_1` and
          `21 delta_0 + (5/2) p_1*p_1`.
        - So in the screened range, B (subspace mixtures) fails only at
          the FM3 endpoint, in the lowest q-degrees.
      - FM-SEC117 (luna_max_pluto; `fm39/sec117_merge_rule_repro.py`): the
        Bernstein step as a merge.
        - KILL: local subgroup merges.  At `(1,2,1,2)`, `j = 1 -> 2`,
          `g_1 = 2 delta_0` has no nontrivial subgroup to merge, while
          `g_2 = delta_0 + (1/2)(delta_1+delta_4)*(delta_1+delta_4)`.
        - KILL: point-by-point origin borrowing.  At `(1^6)`, `j = 0` the
          increment needs 6 units of origin mass against `g_0(0) = 5`.
        - Grouped borrowing into two triangles fixes the increment, but
          the old remainder `g_0 - 3 delta_0` has Walsh minimum -3.
        - Scan of 3,212 Bernstein steps (length `<= 6`): every increment
          is entrywise nonnegative and zero at the origin.
        - Knob: locality of the step; old and new factors must be
          reorganized together.
      - **FM-SEC110 (luna_max_mercury; `fm39/sec110_hACq_d1_repro.py`,
        rerun exactly): H_AC_q is PROVED at distance one, for every list.**
        - Let `n = sum(mu) - 2` and `t` the number of labels 1 in `mu`.
          With `p` the indicator of the single 1's,
          `f_q = ([n]_q!/2) (p*p) + [n]_q! (A_mu(q) - t/2) delta_0`,
          where
          `A_mu(q) = sum_(i<j) q^(mu_(i+1)+...+mu_(j-1)) [mu_i]_q [mu_j]_q`.
        - The factors are independent of `q`, and every coefficient is
          nonnegative: `A_mu`'s coefficients count cross-block endpoint
          pairs at fixed distance, with constant term `m - 1`.
        - Counting reading: the exponent of `q` is the total crossing
          number, `inv(pi)` plus the gap crossed by the extra chord, and
          `[n]_q!` is the inversion enumerator.
        - Exact census: 1,361 ordered lists.  At `(1^8,6)` this reproduces
          FM-SEC102's `a(q), b(q)`.
        - Fixed-factor failures (not failures of H_AC_q): the `q = 0`
          sphere factors at `(1^8, 2)`, where the radius-2 coefficient is
          `-q(1-q)(1-q^2)/30`; and the FM-MECH34 factors at `(1^4, 2)`,
          `n = 2`.
      - FM-CHK50 (luna_max_uranus, own code): the distance-one H_AC_q
        theorem and the `(1^N, n)` closed forms.
        - ACCEPT: the identity, the support argument and the coefficient
          signs (500 random ordered lists; all 48 chord/permutation cases
          of `(1,2,2)`).
        - ACCEPT: the continued fraction and Touchard-Riordan formulas,
          with `M_q = [H_n] H_1^r`.  Raw moments carry an extra `[n]_q!`.
        - REPAIR (conventions only):
          - the per-matching crossing formula `inv(pi) + gap` needs
            right-to-left ranks of the inserted-block legs;
          - the `(1^8,2)` certificate is stated for the
            coefficient-normalized table, and for the raw moment table
            every coefficient is multiplied by `[2]_q! = 1 + q`.
          The sphere-only `b_1` failure holds in both normalizations.
      - FM-SEC112 (luna_max_mars; `fm39/sec112_hACq_spheres_repro.py`,
        rerun exactly): the one-label family `(1^N, n)` under q.
        - Closed forms: `M_q` as a q-Motzkin path sum and continued
          fraction; `C_q(j)` by Touchard-Riordan.
        - The unique sphere coefficients `b_h(q)` fail coefficientwise
          first at `(8,2,3)`, with `b_1 = -q(1-q)^2/30`.
        - H_AC_q holds there by an exact repair: the sphere factors plus
          28 fixed two-point factors `(delta_0 + delta_v)`, `|v| = 6`, with
          weight `(q+q^3)/12`, re-expanded on all 256 masks.
        - The `h >= 2` sphere coefficients are nonnegative for `n <= 12`,
          `d <= n+1`.
      - FM-SEC116 (luna_max_mars; `fm39/sec116_spheres_q_d3_repro.py`,
        rerun exactly): `(1^N, n)` under q, by distance.
        - `d = 2`: spheres alone give coefficientwise nonnegative
          `b_0, b_1, b_2`, for every `n >= 1` (closed forms in
          `S_m = sum [k]_q` and `M_2`).
        - `d = 3`: `b_0, b_1, b_2` are nonnegative for `n >= 3`; `b_3` is
          open.  The case `(8,2,3)` has the two-point repair.
        - Sphere-only census (`N <= 20`, `d <= n + 1`, 69 cases): the only
          failures are at `b_1`, on the line
          `(8,2,3), (11,3,4), (14,4,5), (17,5,6), (20,6,7)`, i.e.
          `d = n + 1`.
        - This is per-`d` progress on one family, so the line is paused
          (full-cone rule).
      - FM-SEC114 (luna_max_uranus; `fm39/sec114_even_repeated_repro.py`,
        rerun exactly): repeated even labels `(2^M)`.
        - A single radial square fails at `M = 5` (all eight sign
          patterns).  B holds there: `f = (2/5) sum 1_(H_ij) +
          (1/5) sum 1_(H_ijk)`.
        - No uniform construction.  The first instance beyond the finite
          certificates is `M = 11`, with profile
          `(1585, 0, 232, 91, 108, 90)`.
      - FM-MECH37 (astra_max_ceres; `fm39/mech37_allocation_repro.py`,
        all assertions pass): channel allocations for the insertion step.
        - KILL: preserving the fusion channel of each noncrossing
          diagram.  At `(1,1,1; h = 1, n = 2)` channelwise CP forces
          `2 d_i <= d <= 1`, so at least half of `(05)(14)(23)` must leave
          its channel 3.
        - The allocation polytope there is exact, with vertices
          `0, (1/2,1/2,0), ...`, each certified by subspaces.
        - KILL: normalized positive Pluecker terminal allocation, even
          after averaging orders.  On `(1^L; 1, L-1)` its channel-`L`
          budget is the harmonic sum `sum_(k<=L) 1/k > 1`.
        - The q-boundary has a clean H_AC_q certificate:
          `g_1(q) = (1+q) Q(p_4)/2 + q(1+q) delta_0`,
          `g_3(q) = (1+q)(1+q+q^2) delta_0`.
        - Unrestricted subset-dependent allocation and H_AC_q survive.
      - **FM-SEC102 (luna_max_mercury): H_AC_q, a unified strengthening.**
        - Conjecture: for every list `lambda`, the q-table
          `f_q(S) = m_q(S) m_q(S^c)` (q-Hermite moments, crossings
          weighted `q`) equals `sum_nu c_nu(q) (p_nu * p_nu)`, with
          factors `p_nu >= 0` independent of `q` and coefficients
          `c_nu in R_(>=0)[q]`.
        - It implies H_AC (at `q = 0`, hence FM3), q-positivity (Q) and
          H_AC at `q = 1`.
        - Screen: all 169 even-total lists with `L <= 7` and labels
          `<= 4` have exact certificates: H_AC at `q = 1`, and
          coefficientwise H_AC_q on 4,346 coefficient vectors (q-degree
          up to 91), each list with one fixed dictionary.  The dictionary
          has subspace-orbit squares, label-block spheres and fusion paths
          `P_(a,b)`.
        - Boundary certificates:
          - `(1^4)`: `f_q = 1_E + (1+q) delta_0`.
          - `(1,5,2,2)`: `f_q = [5]_q! delta_0`.
          - `(1^8,6)`: `f_q = a(q) delta_0 + b(q) (p_7 * p_7)`, with
            `a, b` having nonnegative coefficients and `b` palindromic of
            degree 15.  The main agent re-expanded it exactly on all 256
            elements.
        - The printed code covers the q-table and the certificate check,
          not the dictionary and census driver.
        - FM-CHK49 (luna_max_eris, independent code): ACCEPT the census
          (all 169 lists, 4,346 coefficient vectors) and all four boundary
          certificates.
      - FM-MECH36 (astra_max_ceres;
        `fm39/mech36_insertion_candidates_repro.py`, all assertions
        pass): insertion-stable strengthenings of H_AC.
        - Distance-2 supports add nothing: every H_AC factor of an actual
          table has them, since `f({i}) = 0`.
        - KILL: a common Gram state across external spins.  At
          `mu = (2)`, spins 1, 3, 5 give a non-PSD matrix (witness -1).
        - KILL: translation lifts `C = sum theta tau_s f_j`.  They have a
          proved insertion step, but the hypothesis fails at
          `(1,1,1; h = 1, n = 2)` (budget 1 < 3) and at all three
          boundaries.
        - KILL: scalar channel shares `C_j = gamma_j C`.  At `(1^5, 2)`
          the Walsh bound gives `gamma_1 + gamma_3 <= 29/30 < 1`.
        - A subset-dependent allocation `C_1(S)` (1/3 or 7/9 of `C`)
          repairs that witness, with 16 exact subspace terms.  Bounded
          screen: 167 instances, 10 scalar and 112 translation failures,
          every child Fourier-nonnegative.
        - Open: an insertion-stable subset-dependent allocation.
      - FM-SEC98 (luna_max_pluto; `fm39/sec98_all1_halfplane_repro.py`):
        all-1 sector above the diagonal.
        - Every even-`g` word through length 14 (10,922 words) is
          certified on the whole admissible half-plane.  Claim R covers
          `s <= q`.  Above the diagonal, the tensor Bernstein coefficients
          in `(q, t)`, `s = q + (1-q)t`, are all `>= 0` (2,714,914
          coefficients).
        - Ordered `{1,2}` lists to 6 blocks (59,126 coefficients) and
          `{1,2,3}` lists to 4 blocks also pass.
        - Termwise arguments fail: `ffgg` has a crossing term
          `(q-s)/2 < 0`, `T(fg-gf) = -s(fg-gf)`, and
          `phi(W(ff)W(gg)) = (q-s)/2`.
        - No proof for arbitrary length.
      - FM-CHK43 (luna_max_saturn, own code): ACCEPT FM-SEC66 (39 records,
        295 direct coefficients), FM-SEC69 (Riordan moments `>= 0`) and
        FM-MECH31 Proposition 3.  REPAIR FM-SEC67: the strict cutoffs hold
        for odd `e` only (applied above).
      - FM-CHK44 (luna_max_uranus, own code): ACCEPT FM-SEC71 (the weighted
        certificate covers every real `p >= 0`, `S >= 0`; margin pairs
        rechecked), FM-SEC78 (`r_(p+1)/r_p < 1`, `r_2 = 11/14`), and the
        counting proof that B fails at `(1^8, 6)`.
        - KILL (main agent, `fm39/hAC_sqrt_test.py`): the canonical
          one-term choice `p = WHT^(-1)(sqrt(f_hat))`, i.e. `C^(1/2) >= 0`
          entrywise.  It is negative on 812 of 1,043 label lists (`L <= 10`,
          labels `<= 6`); the worst entry is `-0.47` at
          `(1,1,2,5,5,5,5,5,5,6)`.  Knob: a single symmetric factor, while
          H_AC allows sums.  H_AC certificates must be genuine sums.
        - Open: the general insertion step.  It needs nonnegative,
          pointwise-disjoint factors compatible across the fusion channels
          `CG(h,n)`.  FM-MECH34 is attempting the proof, FM-SEC85
          (luna_max_vesta) the falsification.
      - **FM-SEC77 (luna_max_pluto): q-Gaussian deformation**
        (`fm39/sec77_qgauss_repro.py`, rerun exactly).
        - Replace `U_n` by the monic continuous q-Hermite `H_n(x|q)`
          (`x H_n = H_(n+1) + [n]_q H_(n-1)`), with `x, y` independent
          q-Gaussians.  Then `m_q(S)` sums `q^(crossings)` over the
          matchings with no chord inside a block, and
          `F_q(T) = sum_S (-1)^|S cap T| m_q(S) m_q(S^c)`.  At `q = 0` this
          is FM3.
        - Screens: all 132,788 even-sign profiles (exhaustive to length 6
          and labels 4; `(1^8, 6)`; lists to length 16) have nonnegative
          coefficients in `q`.
        - `q = 1` is proved for every list.  Rotate
          `X, Y = (G +- H)/sqrt 2`; the Hermite addition formula gives
          `F_1(eps) = 2^(L - D/2) sum_(k: (-1)^(n_i-k_i) = eps_i)
          prod C(n_i,k_i) W(k) W(n-k) >= 0`.  This agrees with the
          `q`-polynomials on 264 patterns.
        - The tight word `h_4 hat S_2^2` (labels `(1,5,2,2)`, value 1 at
          `r = 1`) has `F_q / 2 = [5]_q!` exactly (main agent).
        - Stronger than FM3: FM3 is the `q^0` coefficient.
      - **Two-parameter braided family (main agent,
        `fm39/qs_family_screen.py`).**
        - Mixed q_ij-Gaussians with `q_11 = q_22 = q` and `q_12 = s`:
          `F_(q,s)(T) = sum_S eps^S sum_(M_1, M_2) q^(cr(M_1)+cr(M_2))
          s^(cr(M_1,M_2))`, with legs in linear order.  `s = 1` is Pluto's
          line; `(q, s) = (0, 1)` is FM3; `(0, 0)` is free independence.
        - `|s| <= q`, positive by rotation (sketch, to be written out).  In
          the basis `f, g = (e_1 +- e_2)/sqrt 2` the braiding has entries
          `(q+s)/2` (pass) and `(q-s)/2` (flip).  Also
          `e_1^(x)n + eps e_2^(x)n = 2^(-n/2) sum_w (1 + eps (-1)^(#g)) w`.
          So the braided Wick formula is a sum of nonnegative terms when
          `|s| <= q`.
        - `s = 0`, positive for `q >= 0`: colours are constant on the
          components of the chord-and-crossing graph, which gives
          `F_(q,0) = sum_M q^cr(M) prod_C (1 + eps^C)`.
        - Screen: 544 profiles (length `<= 6`, labels `<= 3`, degree `<= 12`)
          on a 41 x 41 grid.  There are no negatives with `q >= 0` and
          `s >= -q`, and negatives exist just below `s = -q`, e.g. at
          `q = 0`, `s = -0.05`.  So the region above the diagonal, which
          contains FM3, is empirically positive, but the rotation does not
          cover it: there the flip weight is negative.  Example:
          `(1,1,1,1)`, `T = {0,1}` gives `F = 4 + 2q - 2s`.
        - FM-SEC90 (luna_max_pluto; `fm39/sec90_braided_repro.py`, rerun
          exactly).
          - Claim R is PROVED: `F_(q,s) >= 0` for `0 <= q <= 1`,
            `|s| <= q`.  The braiding `T(e_i (x) e_j) = q_ij e_j (x) e_i`
            satisfies Yang-Baxter (both sides carry `q_ij q_ik q_jk`) and
            has norm `max(|q|,|s|)`.  In the basis `f, g` its entries
            `(q +- s)/2` are nonnegative, and the signed blocks expand with
            coefficients 0 or `2^(1-n/2)`.
          - Claim Z is PROVED: `F_(q,0) = sum_M q^cr(M) prod_C (1 + eps^C)`
            for `q >= 0`.
          - The boundary is sharp: `(1,1,1,1)`, `T = {0,2}` gives
            `F_(q,s) = 2(q+s)`, which vanishes on `s = -q`.
          - Screens: no negative value in `q >= 0, s >= -q` on 544
            profiles (41 x 41 grid), nor on all 512 masks of
            `(1^8, 2, 2)` (651 grid points each).
          - Further values: `(1,5,2,2)`, `T = {1,2}` gives
            `F = 2 [5]_q!`, independent of `s`.  `(1^8,6)`,
            `T = {0,1,2,3}` gives `F_(0,s) = 24 + 4s - 2s^2 - 8s^3 - 6s^4 -
            4s^5 - 2s^6`, with Bernstein coefficients
            `24, 74/3, 126/5, 126/5, 358/15, 58/3, 6`.
          - No q <-> s duality exists (main agent: the profile sets of six
            lists are not closed under the swap), so Claim Z does not
            transfer to `(0,1)`.
          - Open: positivity above the diagonal, which contains FM3.
        - Bernstein positivity along the free-to-tensor segment (main
          agent, `fm39/qs_bernstein.py`, `fm39/qs_bernstein_fock.py`, and
          the exact Fock-space evaluator `fm39/qs_fock.py`, checked against
          enumeration on 260 values).  For `q = 0`, write
          `F_(0,s) = sum_j b_j C(d,j) s^j (1-s)^(d-j)`, with `d` the degree
          in `s`.  Then every `b_j >= 0` on 542 profiles (length `<= 6`)
          and 2,512 profiles (60 random lists, length 5..8, labels `<= 4`).
          This holds although `F_(0,s)` is not monotone (289 of 542
          profiles increase somewhere).  `b_0` is the free value
          (positive) and `b_d` is FM3.  So this is a strengthening of
          FM3, with the other coefficients interpolating.
        - Reading (main agent): `F_(0,s)` is the limit of averages over
          random bipartite commutation graphs.  Take
          `x = N^(-1/2) sum_a x_a` and `y = N^(-1/2) sum_b y_b`, with
          free families, `x_a` commuting with `y_b` exactly when the edge
          `ab` is present, and edge probability `s`.  The complete graph
          gives FM3 and the empty graph the free value.
        - FM-SEC95 (luna_max_mercury): falsification of the
          strengthenings.
          - `(1^7, 5)`, `T = {0,1}` gives
            `F_(0,s) = 20 + 2s - 2s^3 - 4s^4 - 2s^5`, which is `-120` at
            `s = 2`.  So the half-plane statement needs the bound
            `s <= 1`; the mixed Fock space is defined only for
            `|q|, |s| <= 1`.  Inside the square this profile is positive,
            with Bernstein coefficients `(20, 102/5, 104/5, 21, 20, 14)`.
          - q-positivity (Q): 35,456 more profiles on 41 targeted lists,
            including `(1^N, N-2)`, `(1^N, N-4)` and lists of length 9-12,
            plus the tight atlas.  No negative coefficient.
          - Graph products (G): 1,632 values for selected `N = 3, 4`
            graphs, none negative.
          - Bernstein positivity (BP): no failure; the extended screen is
            incomplete.
        - FM-SEC93 (luna_max_uranus): no proof of Bernstein positivity.
          - Slot reading: `F_(0,s)` is the expectation of the polarization
            `sum_c a_c e_c(z)/C(d,c)` on `d` Bernoulli-`s` slots.  The
            slots are abstract and are not identified with edges or
            crossings.
          - The monotone step fails: `b_1 - b_0 = -2` at `(1^4)`,
            `T = {0,1}`.  The de Casteljau steps only propagate.
          - Screens without failure: fixed-`q` Bernstein at
            `q = 1/4, 1/2, 3/4, 1` (912 profiles each), and bivariate
            Bernstein on the rotation region `|s| <= q` and on the
            half-plane region `-q <= s <= 1` (544 profiles).
        - FM-SEC92 (luna_max_uranus): no transport proof along `s = 1`.
          - A block-product rotation is impossible for `q < 1`, since
            `Cov(U^2, V^2) = (q-1)/2`.
          - A termwise crossing flow fails at `(1,1,2,2)`,
            `T = {0,2}`, which has two summands `-(1+q)` and total
            `2(1+q)^3`.
          - All 190 atlas cases with `F/2 = 1` factor exactly into
            q-integers.
          - Correction to the card: `(1,1,1,1)`, `T = {0,1}` has
            `F_q = 2 + 2q` at `s = 1`.
        - FM-SEC91 (luna_max_jupiter): local involutions fail.
          - The component move is a sign-reversing involution where it
            applies.
          - First-crossing smoothing is not an involution: at
            `(1,1,1,1)` it maps a negative atom to a positive state fixed
            by both rules.
          - At `(1,2,1,2)`, `S = {0,2}`, `D_1 = (03)`,
            `D_2 = (15)(24)`, no smoothing keeps every cluster on one side.
          - Knob: locality of the move.  The exact atlas has 9,586,981
            even-`T` cases (length `<= 8`, labels `<= 4`, ordered), none
            negative.
        - Graph-product family (main agent, `fm39/graph_product_screen.py`).
          - Take free families `x_1..x_N`, `y_1..y_N` and a bipartite
            graph `G`, with `x_a`, `y_b` commuting iff `ab` is in `G`.
            Put `x = N^(-1/2) sum x_a` and `y = N^(-1/2) sum y_b`; both are
            exactly semicircular.
          - The EVEN form `F_G(T)` is FM3 for the complete graph and the
            free value for the empty graph, checked on 84 values.
          - `F_G >= 0` for all 7 graph classes at `N = 2` (2,632 profiles)
            and all 36 classes at `N = 3` (13,536 profiles; lists of
            length `<= 6`, labels `<= 3`, degree `<= 10`).
          - Since `N = 1` with one edge is FM3 itself, this strengthening
            does not by itself reduce FM3.
        - KILL (main agent, `fm39/crossing_set_restriction_kill.py`):
          positivity under an arbitrary set `R` of allowed mutual-crossing
          positions.  `F_R` counts only pairs whose mutual crossings (as
          4-sets of legs) all lie in `R`.  `F_R` is negative for 6,712 of
          688,872 pairs `(profile, R)` (legs `<= 8`).  The witness `-6` at
          `(1,1,1,1,2,2)`, `T = {0,2}`, with an 8-element `R` is exact.  It
          is that profile's exact minimum over all `2^14` sets (FM-CHK47).
          Profiles with more than 14 candidate sets were only sampled, so
          `-6` is not established as the global minimum.  Knob:
          arbitrary crossing sets, which is stronger than any
          colour-based restriction.  So the Bernstein positivity is not
          explained by positivity of every `F_R`.
        - Component model (main agent, `fm39/component_model_check.py`,
          exact on 544 profiles).
          - Group the pairs `(S, M_1, M_2)` by the uncoloured matching
            `M = M_1 u M_2`.  Its components are connected noncrossing
            diagrams, and the colourings are the proper 2-colourings of
            their crossing graph `Gamma`.  So
            `F(T) = sum_(M: Gamma bipartite) prod_(K = (K_1,K_2))
            (eps^(K_1) + eps^(K_2))`.
          - A factor is negative exactly when both sides of a crossing
            cluster are T-odd.
          - Clusters do not cross each other, so this is the free
            moment-cumulant formula: `F = sum_(sigma in NC) prod w_T(K)`,
            where `w_T` is the free cumulant of the commuting variables
            `W_i = U_(n_i)(x) + eps_i U_(n_i)(y)`, with `x, y` commuting
            and independent under the product state (FM-CHK47 repair: for
            freely independent `x, y` the value at `(1^4)` is 4, not 2).  It can be negative
            (`-2` at `(1,1,1,1)`, `T = {0,1}`).
      - FM-SEC87 (luna_max_neptune; `fm39/sec87_even_gf_repro.py`):
        EVEN-form label generating functions at fixed length.
        - For every `L`, `G_T = N_(L,T) / D_[L]` with
          `N_(L,T) = sum_A (-1)^|T cap A^c| K_|A|(x_A) K_|A^c|(x_(A^c))
          D_(A,A^c)`, where `K_m = H_m D_[m]` are the Pluecker numerators
          (`K_4 = 1 - m_1111`; `K_5`, `K_6` explicit).
        - Checked against fusion counts on label boxes for `L <= 6`.
        - Exact positive edge-cone certificates for every even sign class
          at `L = 4, 5`, with 5 to 72 orbit generators.  So FM3 holds for
          every EVEN word of length `<= 5`, all labels (not new: the
          factor axis `<= 6` was proved earlier in the repo).
        - `L = 6` LPs are too large, and the certificates show no rule in
          `L`.  This line is closed.
      - FM-SEC94 (luna_max_jupiter; `fm39/sec94_charging_repro.py`, rerun
        exactly): charging in the component model.
        - Single-smoothing charging fails first at `(1,2,1,2)`,
          `T = {0,1}`: every smoothing of `(03)(15)(24)` leaves an internal
          crossing.
        - Recursive Pluecker straightening reaches positive targets in
          two steps.  Max flow meets the demand on the four boundary
          cases, e.g. demand 38 against capacity 58 at `(2^6)`,
          `T = {4,5}`.
        - Open: the Hall-type capacity condition for the recursive rule.
      - FM-SEC96 (luna_max_vesta; `fm39/sec96_linear_injection_repro.py`,
        rerun exactly): linear injections.
        - The meet identity `U_S cap U_S' = U_(S cap S')` fails at
          `(1^6)` (intersection dimension 1 against 2).
        - Single-flip projections, pair-flip projections and
          `mu_+^* mu_-` all lose rank.  The first loss is at
          `(1,1,1,1,2)`, `T = {0,1}` (rank 3 < 4), where the four negative
          subspaces are linearly dependent inside `H` (kernel
          `(1,-1,-1,1)`).
        - So no map factoring through `H` can be injective.  An injection
          must use the tensor factorization `Inv(S) (x) Inv(S^c)` itself.
      - FM-SEC101 (luna_max_saturn): primal CP search for H_AC.
        - No counterexample or separator.  The boundary certificates are
          re-verified: `(1^8,6)`, `(2^6,4)`, `(1^4,2^3)`, `(1^13,3)`.
        - New: `(2^5,4)` has a clique-orbit certificate, with orbits of
          `{0,3,5,6}` (10) and `{0,3,14}` (60).
        - Of 494 multisets with `L <= 8`, labels `<= 4`: 245 have zero
          targets, and 162 are certified by fusion-path, profile and
          clique dictionaries.  87 are unresolved in those dictionaries,
          which does not mean they fail.  The first is
          `(1,1,2,2,3,3,3,3)`.
        - The random `L = 9, 10` batch and a general rank-`|G|` search
          were not run.
      - FM-SEC99 (luna_max_venus; `fm39/sec99_plucker_grouping_repro.py`,
        rerun exactly).
        - The coset identity `f = sum_M 1_(A_M + W_M)` is proved and
          checked on the boundary atlases (e.g. 20,160 matchings at
          `(1^8, 6)`).
        - The canonical recursive Pluecker family (root plus all
          descendants of first-crossing smoothings) is not H_AC-closed.
          At `(1^6)` the root `(02)(14)(35)` has a 7-member family with
          Fourier coefficient `-1`: 4 of the 10 rooted families there are
          negative, and 7 of 659 at `(1^4, 2^3)`.
        - The whole `(1^6)` table is H_AC (15 subgroups).
        - Knob: grouping granularity.  Groups must be coarser than rooted
          families.
      - FM-SEC100 (luna_max_uranus; `fm39/sec100_level1_census_repro.py`):
        level 1 in the component model.
        - Exact census of `P_T` (i, j on the same side) against `N_T`
          (opposite sides): 5,844,840 ordered profiles (length `<= 9`,
          labels `<= 4`).  `P_T > N_T` in all but 13, where equality holds
          (the degenerate lists); none have `P_T < N_T`.
        - The local one-crossing injection fails at `(1,2,1,2)`,
          `T = {0,1}`, `M = (03)(15)(24)`, as in FM-SEC94.
        - At `|T| = 4` two negative clusters can multiply to a positive
          factor.
        - Open: a nonlocal charging (Hall condition) at level 1.
        - Reproducer repair (main agent): the printed consumer check built
          `h_4` from `U_k(x) U_(4-k)(x)` and printed 25.  With
          `U_(4-k)(y)` it gives `phi_1(h_4 hat S_2^2) = 1`, which agrees
          with `mech28_eval.py`.  The census totals rerun exactly.
      - FM-SEC107 (luna_max_vesta; strategist pass;
        `fm39/sec107_strategist_screens_repro.py`, rerun exactly).
        - Character defect domination:
          `Delta_Q(a,b) = sum_(c in a(x)b) q_(c0) - q_(ab) >= 0` for
          `Q = S_k F_T D_1^(t-2)` and all `a, b >= 1`.  Screened on
          49,954 matrices, none negative.
          - Main-agent reading: `Delta_Q(a,b) = (1/2) E[Q D_a D_b]`, the
            EVEN value of the word with two more minus labels, which has
            `2r + 2` general h.  So this is U-FM3 (the unconstrained cone
            of FM-SEC56) with one extra pair, not a new mechanism.
            `(a,b) = (1,1)` is the consumer value.
        - KILL: a Lee-Yang root condition on the palindromic split counts.
          `(1^4)`, `T = [4]` has roots `+-i` (237 of 3,711 fail).
        - KILL: split-count descent `a_0 >= a_1 >= ...`.  `(1^4)`,
          `T = {0,1}` gives `(3,4,3)` (1,396 fail).
        - No new mechanism emerged.
      - FM-SEC104 (luna_max_vesta; `fm39/sec104_channel_transport_repro.py`,
        rerun exactly).  Equal-label singlet-channel transports on the
        tensor factors, in three variants (sum, signed sum, first pair),
        all lose rank at `(1,1,1,1)`, `T = {0,1}` (rank 1 < 2, kernel
        `(-1,1)`).  General channels need explicit recoupling
        coefficients.  The linear-injection line is paused: its natural
        maps fail on the smallest cases, like the local involutions
        (knob: naturality and locality of the map).
      - FM-SEC103 (luna_max_jupiter;
        `fm39/sec103_terminal_charging_repro.py`, rerun exactly after a
        missing parenthesis was restored): terminal straightening charge.
        - KILL: the charge to fully straightened (noncrossing) targets
          fails.  First at `(1,1,1,1,2)`, `T = {0,2}`: `M = (03)(15)(24)`
          (weight `-2`) straightens to four terms.  One is inadmissible
          (a chord inside the label-2 block) and three have weight 0, so
          the flow is 0.  All 10 failures in the census are permutations
          of this list.
        - At `(1^8, 6)`, `T = {0,1,2,3}` the max flow is 24 against a
          demand of 30.
        - A positive admissible intermediate, `(05)(13)(24)` of weight
          2, exists.  Knob: terminal-only targets.
      - FM-SEC108 (luna_max_uranus;
        `fm39/sec108_fractional_charge_repro.py`, run from the repo root,
        reruns exactly).
        - KILL: canonical fractional straightening charges, in three
          normalizations.  At `(2^5)`, `T = {0,1}` four negative roots
          (weight `-2`, coefficient `c_(M,N) = 2`) all charge
          `N = (09)(12)(34)(56)(78)` (weight 2).  It receives
          `14/3, 16/5, 3888/805`.
        - Also overcharged: `(2^6)`, `T = {4,5}`, and 34, 9, 34 of the 36
          pair-`T` profiles at `(1^8,6)`.
        - Knob: fan-in, i.e. locality of the canonical charge.  Only a
          global flow (Hall) remains on this line.
      - FM-SEC109 (luna_max_jupiter; `fm39/sec109_gamma_star_repro.py`,
        rerun exactly).
        - KILL: charging to intermediate admissible configurations
          (`Gamma*`).  First Hall failure at `(1,1,1,1,2,2)`,
          `T = {1,2}`: the three weight `-2` configurations
          `(02)(17)(34)(56)`, `(07)(13)(24)(56)`, `(07)(15)(26)(34)` have
          the single positive neighbour `(07)(12)(34)(56)`, of weight 4.
          Demand 6 against capacity 4; max flow 12 < 14.
        - Zero-weight relays do not repair it.  `(2^6)` and `(1^8,6)`
          pass (38/38, 30/30).
        - With FM-SEC103 and FM-SEC108 this closes charging by Pluecker
          straightening: terminal, canonical fractional and intermediate
          targets all fail.  Knob: reachability by straightening moves.
      - FM-SEC106 (luna_max_uranus): recursive Pluecker charging at
        level 1.  Max flow meets demand on 1,450 ordered profiles
        (length `<= 5`, labels `<= 3`, all pairs `T`) and on the boundary
        cases.  No uniform Hall argument yet: the overlap of the
        recursively generated target families is uncontrolled.
      - FM-SEC105 (luna_max_venus; `fm39/sec105_groupings_repro.py`,
        rerun after dedenting).
        - Straightening-target and block-partition groupings fail at
          `(1^4)` (Fourier `-1`, or `-1/2` when normalized).
        - The noncrossing-partition grouping factors,
          `G_sigma(S) = prod_(B in sigma) f_B(S cap B)`, and passes the
          screens, but its one-block group is `f` itself.
        - Grouping lines closed.
      - FM-SEC97 (luna_max_mars; `fm39/sec97_spheres_repro.py`, rerun after
        two parenthesis fixes).
        - The sphere expansion for `(1^N, n)` is reproved by Lagrange
          inversion; it has no negative coefficient exactly when
          `d <= n + 1`.
        - Inner range: two-shell factors `p_a + t p_b` certify
          `(1^13,3)`, `(1^16,4)`, `(1^20,2)` exactly.  A fixed ratio
          `t = 1/4` is separated at `(1^28, 8)`; `t = 1/10` repairs it.
        - REPAIR: the printed `(1^30, 2)` coefficients do not re-expand
          exactly (close to the target; transcription).  It is skipped in
          the reproducer.
        - Mixed `d = 2` samples `(1^6,2,4)` and `(2^5,6)` pass.
        - No uniform rule for `d >= n + 2`.
      - FM-SEC85 (luna_max_vesta): no H_AC counterexample.  All 791 tables
        are doubly nonnegative.  Horn 5-cycle placements (374 tables
        exhaustively, 100,000 samples on four larger ones) and a copositive
        `J - 2A` on a 7-vertex triangle-free graph (3,603,600 placements
        for `(1,1,1,1,2)`) give no negative pairing; the tightest
        normalized pairing is 1/2.  Open CP case below the sphere range:
        `(1^13, 3)`.
      - FM-SEC83 (luna_max_jupiter; `fm39/sec83_atlas_repro.py`, rerun
        exactly).
        - The noncrossing-matching model agrees with fusion counts on 3,003
          lists.
        - Atlas: 789,503 even-`T` cases (length `<= 10`, labels `<= 5`),
          none negative.  Every non-parity zero (26) has empty support.
          `F/2 = 1` in 190 cases.
        - A single-cluster crossing switch fails already at `(1,1,1,1)`,
          `T = {0,1}`, `S = {1,2}`: the chords `(1,2)` and `(0,3)` are nested.
        - By hand (main agent): whole-component moves pair the nested
          negatives with one-side positives.  The mutually crossing
          negatives `S = {0,2}, {1,3}` need a crossing-resolution (skein)
          move onto `S = {0,1}, {2,3}`.
      - **FM-MECH31 (astra_max_ceres): B is proved for every
        repeated odd-label list `(q, ..., q)`, all odd `q`, any length**
        (`fm39/mech31_hypB_oddlabels_repro.py`, rerun exactly).
        - General theorem.  Let `Z` be symmetric, `mu_j = E[Z^j]`,
          `d_1..d_L >= 1`, `d(S) = sum_(i in S) d_i`, `D = d([L])` even.
          Put `alpha = min(|Z|,|Z'|)` and `beta = max(|Z|,|Z'|)` for
          independent copies.  Then
          `mu_(d(S)) mu_(D-d(S)) = E[alpha^D] 1_(E_d) + (1/2) sum_(A != 0)
          E[alpha^(D-d(A)) prod_(i in A)(beta^(d_i) - alpha^(d_i))]
          1_(E_d cap C_A)`.
          Here `E_d = {d(S) even}` and `C_A = {bits in A equal}`; all
          coefficients are `>= 0`.  (Expand
          `beta^(d_i) = alpha^(d_i) + (beta^(d_i) - alpha^(d_i))`.)
        - With `Z = U_q(x)`, `q` odd, this gives B for `(q^L)`.
          Insertion of two equal labels propagates positively inside the
          sector.  Consumer content: FM3 for `h_(q-1)^(2r) hat S_q^k` at
          every level (`q` odd).
        - Obstruction: the equality-block family fails at `(1^6, 2)`
          (separating functional, pairing `-15`).  Full B holds there:
          `f = 2*1_0 + (9/2) Avg(1_(E_5)) + (5/2) Avg(1_(K_(2,3)))`, with
          `K_(2,3)` the cycle code of `K_(2,3)`.
        - First open insertion: one label 2 into an all-ones list,
          `sum_(|S| even) y(S) Cat(|S|/2) [Cat(m-|S|/2+1) - Cat(m-|S|/2)]
          >= 0` for every `y` in the dual subspace cone on `F_2^(2m)`.
        - Screens: `(q^L)` for `q = 1,3,5`, `L <= 12`; `(2^k)` and
          `(1,1,2^k)` for `k <= 12`; one label 3 with the rest in `{1,2}`,
          length `<= 8`.  All have certificates.
      - Killed (knob: a restricted subspace family): "B-lite", using
        only the sub-cube subspaces
        `W_T = {X subset T : |X cap Odd| even}`.  It is infeasible on 25
        of 39 random lists (length 4..8, labels `<= 4`;
        `fm39/hypB_lite.py`), so B needs richer subspaces.
    - FM-SEC69 (luna_max_venus), the `{1,2}` sector
      (`hat S_2^k h_1^a`, every level).
      - Proved, with an explicit positive formula:
        `phi_r(hat S_2^k) = (1/2) sum_(p,l) C(k,p) C(2r,2l)
        A_(p,r-l) A_(k-p,l)`, with `A_(p,s) = E[X^(2s) U_2(X)^p] =
        sum_l C(s,l) mu_(p+l) >= 0` and `mu_n = E[U_2^n]` the Riordan
        numbers (noncrossing partitions without singletons), for all
        `r, k`.  This is the parity mechanism made explicit.
      - Even `k` with even `a` is pointwise positive; odd `a` gives 0.
      - Open: odd `k >= 3` with even `a >= 2`, the first mixed-parity
        case.  Separated moments have mixed signs there
        (`(r,a,k) = (1,2,3)`: `7/2, -1/2, -1/2, 7/2`, total 4).
        FM-SEC72.
    - **The H-only level-2 sector as one inequality (main agent,
      `fm39/mp2_kostka_check.py`, with FM-MECH28's evaluator
      `fm39/mech28_eval.py`).**
      - `phi_2(h_(kappa_1) ... h_(kappa_n)) = sum_lambda K(lambda, kappa)
        w(lambda)`, with no bound on `n`, holds with exact equality on all
        137 partitions `kappa` of size 2..10.  Here `w` is the MP_2 weight.
      - On Sp(4), `(x-y)^2 = s_2 - 3 s_(1,1) + 8`, so
        `F_2 = (h_2^perp - 3 e_2^perp + 8) F_1`.  The sector is
        `3 m_(1,1)(E) <= m_(2,0)(E) + 5 m_(0,0)(E)` for
        `E = (x)_i Sym^(kappa_i) W`.  FM51 gives
        `m_(1,1) <= m_(2,0) + m_(0,0)` for every polynomial module.  The
        factor 3 fails for single Schur modules, so a proof must use the
        product structure, as an invariant cone under `(x) Sym^k W`.  This
        is M5 inside one sector with unboundedly many factors.
        FM-SEC60 (luna_max_venus) is on it.
      - Sharpness (`fm39/mp2_ratio.py`, all 271 partitions of size
        `<= 12`).  `max 3 m_(1,1)/(m_(2,0) + 5 m_(0,0)) = 369/386 = 0.956`,
        at `kappa = (3,3,3,1,1,1)`.  No equality case except the 138
        odd-size partitions, where all three multiplicities vanish.
      - The inequality is asymptotically sharp
        (`fm39/mp2_ratio_families.py`).  The ratio tends to 1 along every
        family tested: `(1^n)`, `(2^n)`, `(3^n)`, `(2,1^n)`, `(3^k,1^k)`,
        `(4^k,2^k)`.  For `E = W^(x)2k`,
        `1 - 3m_(1,1)/(m_(2,0)+5m_(0,0)) = 4/(k^2+3k+4)` exactly for
        `k = 1..8`.
        - Reason: `Q_2 = V_(2,0) - 3V_(1,1) + 5V_(0,0)` has virtual
          dimension `10 - 15 + 5 = 0`, since `(x-y)^2` vanishes at the
          identity.  As factors accumulate, the multiplicities approach
          the ratio `1 : 5 : 10`, and `phi_2` is a second-order term near
          the identity.  The same holds at every level
          (`(x-y)^(2r)` vanishes to order `2r` there).
        - Consequence for M5: every consumer value is a
          dimension-zero functional, and FM3 is asymptotically tight as
          factors are added.  A merge-stable hypothesis or invariant cone
          must be exact at leading order, so crude bounds cannot close
          it.  This is consistent with the kills recorded above.
      - FM-SEC60 (luna_max_venus): the sector screen reproduces (915
        partitions of size `<= 16`, none negative, 464 zeros).  The shifted
        cone `{L_(a,b) >= 0}` with
        `L_(a,b) = m_(a+2,b) + 5m_(a,b) - 3m_(a+1,b+1)`, for
        `(a,b) in {(0,0),(1,0),(1,1),(2,0)}`, contains every `Sym^k W` but
        is not preserved: `L_(1,0)(E_(3,3,3,1) (x) W) = -1`, with 70
        failures down to `-329` (knob: translated inequalities, a
        strength the sector does not use).
      - Killed, tentatively (main agent, `fm39/mp2_qkostka.py`; knob:
        coefficientwise `q`-positivity, a strength the sector does not
        use): the Kostka-Foulkes strengthening.  Replacing `K(lambda,kappa)`
        by `K_(lambda,kappa)(q)` (charge statistic) gives polynomials whose
        values at `q = 1` match `phi_2` exactly, e.g. 2 at `(2,1,1)` and 6
        at `(1^4)`.  But they have negative coefficients, e.g.
        `-2 + 3q + q^2` at `kappa = (2,1,1)`.  Tentative because the
        charge code is the main agent's own and has not been checked
        independently.
      - Next route (FM-SEC62).  By Weyl's formula `m_lambda(E_kappa)`
        counts walks in the `C_2` Weyl chamber, one Pieri step per factor.
        So `phi_2(E_kappa)` is a signed walk count, and a sign-reversing
        involution would be exact at leading order.
      - FM-SEC62 (luna_max_venus): `(1^n)` is immediate, since the
        integrand `(x-y)^4 (x+y)^n` is `>= 0` for even `n`.  A rematching
        rule that changes only the last step of a walk fails at length 4
        (Hall deficiency 2 at the prefix `(2,1),(3,1),(4,1),(4,2)`).  The
        even-size zeros through size 16 are exactly
        `kappa_1 > sum_(i>1) kappa_i + 2`.  General `Sym^k` factors need
        the actual Pieri transitions.
      - Killed (main agent, `fm39/mp2_local_kostka.py`; knob: local
        matching, stronger than the global sum): the shape-by-shape route.
        Each negative shape `lambda = (k+b+2, k+b+1, k+1, k)` has the
        positive neighbours `lambda + e_1 - e_3`, `lambda + e_2 - e_4` and
        `lambda + e_2 - e_3` (weights 1, 1, 2).  Each positive shape of
        weight 1 or 2 is the neighbour of exactly one negative shape.  But
        the local inequality
        `2K(lambda,kappa) <= K(lambda+e_1-e_3,kappa) +
        K(lambda+e_2-e_4,kappa) + 2K(lambda+e_2-e_3,kappa)` fails in 217
        of 461 cases with size `<= 14` (first at `lambda = (3,2,1)`,
        `kappa = (2,1,1,1,1)`: `16 > 12`).  So the weight-5 shapes and the
        global structure must take part.
    - FM-SEC52 (luna_max_neptune), level `r = 1` on the whole cone:
      exact screen of 676,368 words (2..6 `hat S` with labels `2..8`,
      at most two `h`, `a <= 10`), none negative, least positive value 1
      (e.g. `h_4 hat S_2^2`).  The merge correction can be negative
      (`phi_1(X_22 h_1^2) = -2` while `phi_1(hat S_2^2 h_1^2) = 2`).  The
      first open case, two `hat S`, is
      `A_p((q-1,1); E) <= A_p((q); E) + 2 A_p((q-2); E)` in FM51's
      notation, for `E = Sym^u W (x) Sym^v W (x) W^(x)a`.  Assigned as
      FM-SEC61.
    - FM-SEC61 (luna_max_neptune): coverage of the level-1 two-`hat S`
      case.
      - S2-2 and S2-3 prove it for `Sym^u (x) Sym^v` and
        `Sym^u (x) Sym^v (x) W`, FM42 for `W^(x)a`, Lemma BP when a part is
        `>= P+Q-1`, and fixed pairs `(P,Q)` for all `kappa`.
      - It fails for single Schur modules: `S_(2,1,1)`, `P = Q = 2`, gives
        `Delta = -1`.
      - Screens: 43,316 consumer cases (`2 <= Q <= P <= 8`,
        `u, v <= 12`, `a <= 16`) and 19,305 general products, none
        negative, least positive value 1.
      - Residual: `a + [u>0] + [v>0] >= 4` and
        `2 <= max(u,v) <= P+Q-2`, outside the fixed pairs.
      - FM-SEC66 extends the S2-2/S2-3 generating-function certificates
        (positive numerator decompositions over the Sp(4) invariant rings)
        to an unbounded suffix.
      - FM-SEC66 (luna_max_neptune) result
        (`fm39/sec66_two_hat_a2_repro.py`, rerun exactly).  The suffix
        `a = 2` case holds for all labels: the numerator
        `N_2 = [z_1 z_2] N_4` over `D_4 = prod (1 - x_i x_j)` has a 39-term
        positive orbit decomposition, checked against direct evaluation in
        295 coefficients.  So the level-1 two-`hat S` case holds for all
        labels at `a = 0, 1, 2`.  For each fixed core, LS covers
        `a + 6 >= Theta = u(u+4) + v(v+4) + (5/3)(P(P+2) + Q(Q+2))`, which
        leaves a finite window; four cores are certified completely.
        Open: `3 <= a < Theta - 6`, uniformly in the labels.
      - FM-SEC79 (luna_max_neptune; `fm39/sec86_two_hat_a6_repro.py`,
        rerun exactly after a bracket typo fix).  The level-1 two-`hat S`
        case holds for all labels at `a = 3, 4, 5, 6`.  The label
        generating function is `N_a / D_4`, and each numerator has a
        positive edge-orbit decomposition, re-expanded exactly: 43, 83,
        116 and 183 terms at degrees 9 to 12.  At `a = 3, 4` the numerator
        was checked against direct Catalan values on 1,001 label tuples
        each.  Open: `7 <= a < Theta - 6`.  At `P = Q = n`, `u = v = 2` the
        LS threshold is `Theta - 6 = 18 + (10/3) n (n+2)`, so no fixed
        suffix cap suffices.  The extraction recurrence for `J_m` has a
        subtraction and gives no positive map `N_a -> N_(a+1)`.
    - **Theorem LS-hat (FM-SEC63, luna_max_jupiter; ACCEPTED by FM-CHK42,
      luna_max_vesta, own code).**  For a fixed core `H = prod h_(kappa_l)`, maximum
      `hat S` label `P`, and `k >= K(r, kappa, P)` factors `hat S`
      (explicit, very large), `phi_r(H prod hat S_(p_i) h_1^a) >= 0` for
      every `a`.  Proof: away from the four corners the product of `hat S`
      decays like `exp(-c_P k tau)`; near the corners every factor keeps
      its corner sign; a positive patch near `(2,2)` gives the lower
      bound; LS covers large `a`.  Consequence: with the labels capped, each
      level has only finitely many words left (LS bounds `a`, LS-hat
      bounds `k`, `m <= 2r`).  Without a cap the residual is infinite,
      e.g. `phi_1(hat S_p hat S_(p+1) h_1^(2p+1))`.
    - FM-SEC78 (luna_max_jupiter): that family is positive for every
      `p >= 2`.  Merge `hat S_p hat S_(p+1)` into
      `sum_j hat S_(2p+1-2j) + X_(p,p+1)`.  The one-`hat S` terms are
      `>= 0` by Theorem OL, and the `hat S_1` term
      `L_p = 12 C(2p+2,p+1)^2 (2p+3)/((p+4)(p+3)^2(p+2)^2)` dominates the
      explicit cross term.  The next residuals are
      `phi_1(hat S_p hat S_(p+2) h_1^(2p))` and
      `phi_1(h_t^2 hat S_(t+2) h_1^(t+2))`.  Family certificates of this
      kind are not pursued further; the full cone needs a mechanism.
    - FM-SEC31 (luna_max_pluto): the mixed term expands exactly through
      lower-label values, `C = A - B`, so the inductive step is
      `A <= R + B`.  Every term there has fewer labels, but positivity of
      the lower words does not compare them.  A Cauchy balance and an LD
      envelope were screened and are not enough.  The base would need
      two-label results for the same kernel rows, not only Theorem OL.
  - **Three labels with `hat S` (M3; FM-SEC23, luna_max_neptune;
    verified by the main agent, `fm39/sp3_repro.py`,
    `fm39/sp3_band_check.py`).**
    - *Patterns.*  At level `r`, with `e = 2r - m` (`m` = number of general
      `h` factors), the three patterns are
      `h_u h_v hat S_p` (`e = 2r-2`), `h_u hat S_p hat S_q` (`e = 2r-1`) and
      `hat S_p hat S_q hat S_s` (`e = 2r`).  All are consumed at every
      `r >= 1`.
    - *Closed forms.*  Exact split formulas in `W_c`, checked in 5,180
      cases against Catalan moments.
    - *Support branch.*  Sort the labels `X >= Y >= Z` and put
      `alpha = (N+X-Y-Z)/2`, `gamma = (N+X+Y-Z)/2`.  If `gamma >= N+1`,
      the value is `sum_(d in CG(X,Y)) [D_j - D_i + sigma W](d, Z)`.  So
      (E) implies positivity there, the analogue of Theorem G0E.
      - The sign depends on the roles of the labels, not only on the
        pattern (repair from FM-CHK40).  Tag `h` labels `-` and `hat S`
        labels `+`, and sort by decreasing label, `-` first on ties.
        - `h h hat S`: `sigma = +1` if the top two tags are both `-`,
          else `-1`.
        - `h hat S hat S`: `sigma = -1` if the bottom tag is `-`, else
          `+1`.
        - `hat S hat S hat S`: `sigma = +1`.
      - Witness: at `r = 1`, `a = 3`, `h_3 h_3 hat S_3` and
        `h_3 h_2 hat S_4` share sorted labels `(4,4,3)` but need opposite
        signs (values 15 and 13).  If also `0 <= alpha <= N - Z`, it
      telescopes to `S_Z(alpha) + sigma T_Z(alpha)`.
    - *Unconditional band* (G0B method).  The band is
      `(a-e)^2 <= 2(N+1)` when `Z >= 3`, and `(a-e)^2 <= N+1` when
      `Z = 2`.  In it, `S > |T|`.
    - FM-CHK40 (luna_max_saturn, own exact code):
      - ACCEPT: the closed forms (5,180 comparisons), the telescoped form
        (1,307 windows) and the band (247,711 wide-band plus 9,818
        `Z = 2` windows).
      - REPAIR (applied above): the role-ranked sign rule.
    - *Checks.*
      - Agent screen: 48,552 words (`r <= 6`, `a <= 16`, labels
        `2..8`), none negative.
      - Main agent: 3,115 band support words across the three patterns,
        none negative.
    - For `sigma = -1` the support value is exactly a W window, so W's
      criteria apply.  For `sigma = +1` it is `sum D + T`, which is
      automatic when `T >= 0`.  The remainder outside the support branch
      is open.
    - FM-SEC34 (luna_max_neptune): exact merge forms off the support
      branch.
      - Census (`r <= 6`, `a <= 16`, labels `2..8`): 19,643 words, all
        nonnegative.  Some grouping has a favourable correction sign in
        about 75% of them.
      - Correction-resistant words (still positive) include
        `h_2^2 hat S_4 h_1^4`, `h_3 hat S_3^2 h_1^3` and `hat S_2^3 h_1^2`
        at `r = 1`.  These lie outside an (E)-only proof.
    - FM-SEC46 (luna_max_neptune): four factors beyond G0E4.
      - For each pairing `pi`, `phi = M_pi - C_pi`, where `M_pi` is a sum
        of (E) values (`d + d' <= N`) and one-sided OL terms.  So
        `phi >= 0` whenever some `C_pi <= 0` and the consumed (E) calls
        are proved.  This branch is not exhaustive: all three `C_pi > 0`
        at `(r, a, L) = (2,3,(8,6,6,3))`, `(3,1,(7,6,3,3))`,
        `(4,3,(9,9,8,3))` (shifted labels; `phi = 54, 18, 263`).
      - A two-label telescope writes each correction as a signed wedge
        sum (4,800 exact checks).  Bounding it by absolute values loses
        the cancellations.
      - Its first unresolved inequality, `C_pi <= M_pi` for some
        pairing, is `phi >= 0` itself.  No reduction.
  - **GFM3: FM3 with extra nonnegative kernel factors (main agent,
    `fm39/recip_test.py`, `fm39/recip_bounds.py`, `fm39/recip_quartic.py`,
    `fm39/kernel_identity.py`, `fm39/general_row_words.py`).**
    - *Identity.*  For `P = (1+z)^a (1-z)^e prod_i (1 + t_i z + z^2)` the
      quadratic kernel formula `W_c(p,q) = B_i c_j - c_i B_j` equals
      `E[K U_p(x) U_q(y)]` with
      `K = (x+y)^a (x-y)^e prod_i Q_(t_i)`,
      `Q_t = x^2 + t x y + y^2 + t^2 - 4`,
      because `P(u/v) P(uv) = u^N K(u+1/u, v+1/v)`.  Checked in 8,692
      cases, any real `t`.
    - `Q_t = (x + t y/2)^2 + (t^2/4 - 1)(4 - y^2)` is `>= 0` on
      `[-2,2]^2` iff `|t| >= 2`.  `Q_2 = (x+y)^2`, `Q_(-2) = (x-y)^2`, and
      `Q_t / t^2 -> 1` as `t -> infinity`.
    - *Screens* (none failed):
      - the three-factor closed form, all `gamma` (including
        `gamma <= N`), on anti-reciprocal real-rooted rows (`e` odd,
        `|t_i| >= 2`): 342,495 tests;
      - OL on these rows: 11,651 tests;
      - (E) on reciprocal rows: 14,250 tests;
      - the full split functional on 20,000 random consumer words
        (2..6 general labels, 0..2 `hat S`, 1..3 factors `Q_t`,
        `|t| >= 2`).
      - a stress screen of 3,000 words with clustered small labels
        (`base, base+1`, `base` in 3..5), 3..8 general factors, the
        minimal `e` (the FM2 boundary), 0..3 `hat S` and 1..4 factors
        `Q_t` (`fm39/gfm3_hard.py`): no negative value and no zero.
    - *Controls* (failures found):
      - `|t| < 2`: 102 of 3,000 words, and 1,291 three-factor failures;
      - reciprocal complex root quadruples (a positive weight, but not
        real-rooted): 17 failures;
      - non-reciprocal real-rooted rows: 9 failures in 49,513.
    - So FM3 appears to be the vertex case `t_i in {+-2}` of a
      positivity statement on the whole family `Q_t`, `|t| >= 2`.  In
      `s = 1/t in [-1/2, 1/2]`, `Q_t / t^2 = 1 + s xy + s^2 (x^2 + y^2 - 4)`.
      The family deforms `(x+y)^2` (`s = 1/2`) through `1` (`s = 0`) to
      `(x-y)^2` (`s = -1/2`).
    - *What the positivity uses* (`fm39/window_ulc.py`,
      `fm39/logconcave_test.py`):
      - Positive rows: Conjecture W held on 600 ultra-log-concave rows
        `C(n,k) u_k`, `u` log-concave, most of them not real-rooted.  Plain
        log-concave rows fail (541 of 600).  So for positive rows Newton's
        inequalities appear to suffice.
      - Mixed signs: rows satisfying Newton's inequalities with signs fail
        in 22 of 600 (random signs) and 7 of 600 (one sign change).  So the
        mixed-sign case uses real-rootedness beyond Newton.
      - The unit-circle weight of the row is `K(x,2)`, since
        `P(u)^2 = u^N K(u+1/u, 2)`.  A factor `Q_t` contributes `(x+t)^2`.
        For complex `t = al + i be`, the factor is `((x+al)^2 + be^2)`,
        which is log-concave on `[-2,2]` iff `|al| - 2 >= be`.  Three-factor
        screen: log-concave weight, 107,862 tests, 0 failures;
        non-log-concave weight, 3 failures in 194,366.  The failures track
        log-concavity of the weight in `x = 2 cos theta`.
    - *Geometric form of Conjecture W* (`fm39/polygon_lemma.py`).  Put
      `g_k = (c_k, c_(k-1))` in the plane.  Then `D_k = g_k ^ g_(k+1)`,
      `T(x, x+C) = g_x ^ g_(x+C+1)`, and `P_C(x)` is twice the signed area of
      the closed polygon `g_x, g_(x+1), ..., g_(x+C+1)`.
      - `D_k >= 0` says the curve turns counterclockwise about `O`.
      - `P_1(k-1) = D_(k-1) + D_k - T(k-1,k) >= 0` says it turns left at
        `g_k` (local convexity).
      - Planar screen: for left-turning counterclockwise polygonal arcs,
        closing chords gave no negative area in 263,887 sub-arcs sweeping
        at most `2 pi` about `O`, and 2 negative areas in 111,336 sweeping
        more.
      - So W splits into: Newton (`D >= 0`); the four-coefficient case
        `C = 1`, which Rolle's theorem reduces to real-rooted cubics; a
        convexity lemma for sweeps `<= 2 pi`; and a separate argument for
        longer sweeps.  Longer sweeps occur: near `a = e` the rows turn by
        about `pi/2` per step.
    - **Theorem WS (Conjecture W for short sweeps; main agent,
      `fm39/window_sweep.py`; repaired after FM-CHK34).**  Let `c` be the
      row of a real polynomial of degree `N` with only real roots and
      `c_0 c_N != 0`.  Let `0 <= x` and `x + C + 1 <= N + 1`.  If the phase
      polygon `g_x, ..., g_(x+C+1)` sweeps a total angle
      `Theta = sum_k arg(g_k -> g_(k+1)) <= 2 pi` about `O`, then
      `P_C(x) >= 0`.
      - *Nondegeneracy.*
        - Adjacent zero coefficients cannot occur.  A repeated zero of a
          derivative of a real-rooted polynomial is a zero of the
          polynomial itself, so iterating back would force `c_0 = 0`.
          Hence every `g_k != 0` on `[0, N+1]`.
        - The strict Newton factor `(1+1/k)(1+1/(N-k)) > 1` gives
          `D_k(c) > 0` for `0 <= k <= N`.  The same holds for
          `f = (1-z)c`, which is real-rooted with `f_0 = c_0 != 0`.
        - So every step angle lies in `(0, pi)`, `Theta > 0`, every
          vertex is a strict left turn, and no edge lies on a ray through
          `O`.
        - Without these hypotheses (only `D >= 0`) zero vectors can
          occur, e.g. `c = (1,0,0,1)`, and the angles are undefined
          there.
      - *Proof.*  Each step turns counterclockwise about `O`, since
        `g_k ^ g_(k+1) = D_k(c) >= 0`.  Each vertex is a left turn, since
        `(g_(k+1) - g_k) ^ (g_(k+2) - g_(k+1)) = f_(k+1)^2 - f_k f_(k+2)
        = D_(k+1)(f) >= 0`.
        - If `Theta < pi`, the polygon `O, g_x, ..., g_(x+C+1)` is
          star-shaped from `O`, hence simple.  It turns left at every
          vertex: at `g_x` and at the last vertex by `D >= 0`, and at `O`
          because `T(x, x+C) = g_x ^ g_(x+C+1) > 0`.  So it is convex.  The
          chord `g_x g_(x+C+1)` splits it into the triangle
          `O, g_x, g_(x+C+1)` and the window polygon, both positively
          oriented, so `P_C(x) >= 0`.
        - If `pi <= Theta <= 2 pi`, then
          `T(x, x+C) = |g_x| |g_(x+C+1)| sin Theta <= 0`, and
          `P_C(x) = sum D_k - T >= 0`.  QED.
      - *Consequence.*  On the G0 branch, Theorem WS proves
        `phi_r(h_u h_v h_w h_1^a) >= 0` whenever the window at `alpha`
        sweeps at most `2 pi`.  This is uniform in `r` and `a`, needs no
        root bounds, and holds for every real-rooted row.
      - *Coverage.*  79,442 of the 533,190 G0 words with `r <= 8`, `a < 40`
        (14.9%).  The phase curve usually spirals several times around
        `O`.
      - *Obstruction to the naive extension* (knob: arbitrary
        left-turning ccw curves, no consumer use).  Splitting at an
        interior vertex gives `P = P_first + P_second + 2 Area(g_x, g_y,
        g_end)`.  The triangle needs a vertex strictly right of the chord,
        on the side away from `O`, and a contracting spiral can have none.
        The random-curve failures are exactly of this kind.  Long sweeps
        need a quantitative input from the row, for example Newton for
        `(1 + t z) P` for every real `t`, or the recurrence of the base
        rows.
    - **Long sweeps: pass criterion (main agent, `fm39/window_passes.py`,
      `fm39/window_cover.py`).**  Suppose `T > 0` and `Theta > 2 pi`, and
      put `I = [theta_x, theta_x + beta]`, `beta = Theta mod 2 pi`.  Every
      direction in `I` is crossed by the first pass (starting at `g_x`) and
      by the last pass (ending at `g_end`).
      - On each edge `u = 1/rho(theta)` is exactly sinusoidal, and left
        turns give `u'' + u >= 0`.  On an interval shorter than `pi` the
        maximum principle then puts the chord `g_x g_end` inside the
        first-pass fan if `rho_first(theta_end) >= |g_end|`.  It puts the
        chord inside the last-pass fan if `rho_last(theta_x) >= |g_x|`.  In
        either case `T <= sum D_k`.
      - More generally it suffices that for some `theta* in I` both passes
        are outside the chord.  Put `Z` on the nearer pass at `theta*`, so
        both pass radii there are `>= |Z|`.  The same chord comparison on
        the two shorter intervals contains the triangles `O, g_x, Z` and
        `O, Z, g_end`, whose signed areas add to that of `O, g_x, g_end`.
      - The radial-graph form `u = 1/rho` needs edges off the rays through
        `O`, i.e. `D_k > 0`.  This holds under the nondegeneracy of
        Theorem WS.
      - FM-CHK34 (luna_max_venus, own code):
        - ACCEPT: the phase-polygon and left-turn identities, the twist,
          and the coverage figures (floating-point screen).
        - REPAIR (applied above): the `Theta = 0` case, the zero-vector
          and collinear conventions, and the radial-graph assumption.
        - Structure of the uncovered words:
          - 94 distinct windows, all at `r = 4..8`, with `T > 0` and
            sweep `> 2 pi`;
          - mostly near-equal top labels: `u - v = 0` (832 words),
            `1` (292), `2` (135);
          - window lengths `C` from 3 to 13.
      - The twist `c_k -> (-1)^k c_k` keeps every `D_k` and keeps
        real-rootedness, and sends each step angle `theta_k` to
        `pi - theta_k`.  It keeps `P_C` for even `C`.  For odd `C` it gives
        `P_C = sum D + T~`.
    - *Coverage of the 533,190 G0 words:*
      - `T <= 0` (trivial): 56.36%;
      - pass criterion: 41.55%;
      - convex short sweep: 1.78%;
      - twisted versions: 0.07%;
      - uncovered: 1,259 words (0.24%), for example
        `(r, a) = (4, 14)`, `alpha = 7`, `C = 5`.  These are
        near-symmetric windows, where the first and last passes nearly
        coincide and the middle full turn carries the area.
    - These are window-by-window criteria, checked in floating point.  A
      uniform proof needs analytic conditions on the row under which one
      of them holds.
    - *Newton certificates (main agent, `fm39/window_lp.py`).*
      - Exact identities:
        `P_1(x) = D_(x+1)[(1-z)P]` and `P_2(x) = D_(x+2)[(1-z^2)P]`.  So
        Conjecture W holds uniformly for `C <= 2`, for every real-rooted
        `P`.  `C = 2` is `w = 1`, where `h_w = h_1` can also be absorbed
        into the suffix.
      - For `C = 3, 4` the LP finds no nonnegative combination of Turan
        determinants `D_k[R P]` with `R` in:
        - `(1-z)^j (1+z)^l`, `j + l <= 6`;
        - these times `(1 + t z)`, for 10 values of `t`.
        Terms were allowed outside the window provided they cancel
        (up to 1,834 columns); the LP was infeasible.
      - So the genuine three-factor windows (`C >= 3`) need more of
        real-rootedness than Newton's inequalities for real-rooted
        multiples: Hermite/Bezout positivity, or the phase-polygon
        geometry.
      - Also killed (knob: an auxiliary induction device, no consumer
        use): `G_C(s) = sum_x P_C(x) s^x` is not real-rooted
        (`fm39/window_gf.py`).
      - The two-label inequality (E) has the same status
        (`fm39/e_lp.py`): with reciprocity imposed, 130 of 135 cases
        (`N <= 12`) have no such certificate.  Multipliers up to degree 4,
        including `(1 + t z)` factors, were tried.
    - **FM-MECH23 (astra_max_ceres; reproducer rerun by the main agent,
      `fm39/mech23_repro.py`).  No movement on W; three exact results.**
      - Identity, for every `P`:
        `P_C(x) = sum_(i=0..floor((C-1)/2)) D_(x+C)[(z^i - z^(C-i)) P]`.
        For `C >= 3` some multipliers have roots on the unit circle.
      - For `C = 3` the only covariance representing `P_3` as
        `sum lambda D[R P]` is `G = I - J`.  It forces
        `R = (1-z)(a + (a+b)z + a z^2)`, whose averaged discriminant is
        `-2 < 0`.  So no fixed positive average of real-rooted multipliers
        works.
      - Degree induction:
        `P_C(x; (1 + rho z)P) = P_C(x) + rho l_C(x) + rho^2 P_C(x-1)`, with
        `l_C(x) = sum_(k=x..x+C) T(k-1,k) - c_x c_(x+C-1) + c_(x-2) c_(x+C+1)`.
        So it closes iff `l_C(x)^2 <= 4 P_C(x) P_C(x-1)` (31,880 exact
        checks, unproved).  Local Newton and cubic inequalities do not
        suffice: the row `(4,-3,1,1,-2,-4)`, which has only one real root,
        satisfies them and has `P_3(1) = -2`.
      - Q_t deformation `F(s) = A + B s + C s^2`: endpoint positivity
        controls the interior unless `C > 0` and `|B| < C`, where
        `4AC >= B^2` is needed.  The shortcut `|B| <= 4A` fails on
        `phi_3(h_21 h_4 h_2 h_1^15)`, with `(A,B,C) = (1843, 9097, 16813)`.
    - **W reduces to a crude energy bound on long sweeps (main agent,
      `fm39/uncovered_slack.py`, `fm39/long_energy2.py`,
      `fm39/long_energy3.py`).**
      - `C <= 2` is proved by the certificates, and sweeps `<= 2 pi` by
        Theorem WS.
      - What remains is `sum D_k >= T` on windows with `C >= 3`, `T > 0`
        and sweep `> 2 pi`.  These carry large slack: the least `sum D/T`
        is 1.48 over all 39,935 such windows (grid `r <= 8`, `a < 40`, any
        `C >= 2`).
      - Sufficient condition LE: `sum_(k=x..x+C) D_k >= |g_x| |g_(x+C+1)|`,
        which implies `sum D >= T`.  It held on the grid (39,746 windows,
        `C >= 3`, least ratio 1.55) but is FALSE beyond it (FM-MECH24).
        Examples: `a = 0`, `e = 21`, `x = 9`, `C = 4`.  The family
        `a = 0`, `e = 2n+1`, `x = n-1`, `C = 4` has
        `S/(|g_x||g_end|) ~ 5/n`.  Also `(a,e,x,C) = (4087, 9, 2003, 90)`.
        Kill knob: LE drops the endpoint angle, which is strength the
        consumer does not use; the target `sum D >= T` still holds there.
      - `|g_k|^2` is not unimodal on 209 of 280 base rows, so a radius
        monotonicity argument is not available.
    - **Theorem G0B (the G0 branch in a growing band; FM-MECH24,
      astra_max_ceres; verified by the main agent, `fm39/g0b_repro.py`,
      `fm39/g0b_check.py`).**  If `(a - 2r + 2)^2 <= 8r - 9`, i.e.
      `(a - e - 1)^2 <= 4e + 3`, then `phi_r(h_u h_v h_w h_1^a) >= 0` on the
      whole G0 branch, uniformly in the three labels.  This includes
      root-crossing windows and both parities of `C`.
      - *Energy.*  `D` is symmetric and unimodal (Theorem OL and
        anti-reciprocity; `D_n = D_(n+1)` when `N = 2n`).  So for
        `S = sum_(k=x..y) D_k`, `y = x+C`:
        `S >= max(D_x, D_y) + C min(D_x, D_y) >= 2 sqrt(C D_x D_y)`.
      - *Chord.*  Put `d = a - e` and
        `rho_k = |d| / (2 sqrt((k+1)(N-k+1)))`.  From the recurrence,
        `D_x = p^2 - d p q/(x+1) + b q^2 >= (1 - rho_x)(p^2 + b q^2)`, with
        `p = c_x`, `q = c_(x-1)`, `b = (N-x+1)/(x+1)`.  Likewise
        `D_y >= (1 - rho_y)(v^2 + b' u^2)`, with `v = c_y`, `u = c_(y+1)`,
        `b' = (y+1)/(N-y+1)`.  Since `b b' >= 1`, Cauchy--Schwarz gives
        `T^2 = (p v - q u)^2 <= D_x D_y / ((1 - rho_x)(1 - rho_y))`.
      - *Window region.*  If `rho_x, rho_y < 1` and
        `4C(1 - rho_x)(1 - rho_y) >= 1`, then `S >= |T|`.
      - *The band.*  `d^2 <= 2(N+1)` gives `rho <= 1/sqrt 2`, hence
        `T^2 <= (6 + 4 sqrt 2) D_x D_y < 12 D_x D_y <= 4C D_x D_y <= S^2`
        for `C >= 3`.  Windows with `C <= 2` are covered by the
        certificates.
      - *Checks.*  The agent verifier reruns (57,014 band windows, 27,276
        of them long-sweep with `T > 0`).  The main agent's
        definition-level evaluator found 393,947 band G0 words
        (`r <= 13`, `a < 60`), none negative.  The main agent also
        re-derived each step.
      - *Open:* windows outside the region, in particular `a` far from
        `e` (`rho` near 1).
    - **Centered windows, a metric chord bound, and the equal-label sector
      (FM-MECH25, astra_max_ceres; verified by the main agent,
      `fm39/mech25_repro.py`, `fm39/eqlabel_check.py`).**
      - *Sign-block bound.*  If the folded row (recurrence parameter
        `delta = |a-e|`) has a block of `L` same-sign coefficients bounded
        by opposite or zero signs, then
        `delta^2 <= (N+1)^2 (1 - 4/(L+1)^2)`.  This follows from the
        largest eigenvalue `(N+1) cos(pi/(L+1))` of the normalized
        tridiagonal recurrence.
      - *Theorem (centered windows).*  For every base row with `e` odd,
        `2x + C = N` implies `P_C(x) >= 0`.  Wide windows
        (`C^2 >= N+4`) follow from a quadratic in `t = c_x/c_(x-1)`;
        narrow windows from Theorem WS or the sign-block bound.
        Equivalently, the equal-largest-label sector is settled:
        `phi_r(h_u^2 h_w h_1^a) >= 0` whenever `2u >= a + 2r + w - 2`, for
        every `r >= 2`, `a >= 0`, `u >= w >= 1`.
      - *Metric chord bound.*  With the quadratic form
        `H(p,q) = ((e+1)(p+q)^2 + (a+1)(p-q)^2)/2` of determinant
        `V = (a+1)(e+1)`, one has
        `T^2 <= K D_x D_y / ((2 sqrt V - R_+)(2 sqrt V - R_-))`, where
        `h = 2x+C-N`, `R_+ = C+h`, `R_- = |C-h|` and
        `K = (N+2-C)^2 - h^2`.  With the energy bound `S^2 >= eta D_x D_y`
        this proves W whenever `eta (2 sqrt V - R_+)(2 sqrt V - R_-) >= K`,
        a region that includes arbitrarily large `|a-e|`.
      - *Decomposition.*  A center-crossing window splits as
        `2 P_C(x) = Z(l) + Z(x) - Delta(l,x)`, where
        `Z(t) = P_(N-2t)(t)` is a centered window,
        `Delta = (c_l - c_x)^2 - (c_(l-1) - c_(x-1))^2` and `l = x-h`.
      - *Checks.*  The agent script reruns: 83,300 centered windows and
        1,375,467 metric checks.  The main agent's evaluator found
        223,715 equal-largest-label G0 words, none negative.
      - FM-CHK38 (luna_max_saturn, own exact code): ACCEPT on all five
        items:
        - the sign-block bound (4,456 blocks);
        - the centered-window theorem, every branch;
        - the equal-label translation (717 identities);
        - the metric chord bound (128,018 windows; the region reaches
          arbitrarily large `|a-e|`, e.g. `a = 0`, `N = 5001`);
        - the centered decomposition (16,880 windows).
    - **W on the grid is now complete** (main agent,
      `fm39/w_residual4.py`).  With the centered theorem and the metric
      criterion added, all 233,865 consumer windows (`r <= 10`, `a < 60`)
      are covered by proved per-window criteria:
      - `T <= 0`: 53.4%;
      - metric criterion: 26.7%;
      - LD plus energy: 9.8%;
      - plateau: 5.7%;
      - `e = 1`: 3.5%;
      - centered: 0.8%;
      - G0E rows: 0.1%;
      - Theorem WS: 0.03%;
      - residual: 0.
      So the G0 branch is proved on this grid.  Open: that the union of
      criteria is exhaustive for all `r` and `a` (sampled scan
      `fm39/w_residual_sample.py`; FM-SEC27).
      - Sampled scan beyond the grid: 2,000 random rows (`r <= 40`,
        `a <= 400`), 34,598,270 consumer windows.  Exactly one window is
        uncovered: `(r,a,x,C) = (35,5,46,4)`, with `u = 3.93` just below
        `C = 4`.  It has `S/T = 7.5` and (ii)-ratio 1.55.  `D` decreases
        nearly geometrically there, and the crude bound `max + C min`
        misses the LD bound by about 1%.  A sharper energy bound on the
        decreasing side would cover it (FM-SEC27).
      - Systematic scans (main agent, `fm39/w_residual_rows.py`, `..._rows2.py`,
        `..._rows3.py`).
        - Without the twist, rows with small `a` (3..5) and large `e`
          (`r >= 24`) leave 707 uncovered windows.  These rows alternate,
          so the phase curve turns nearly `pi` per step.
        - Adding Theorem WS on the twisted row (even `C`), the odd-`C`
          twist test `T~ >= 0`, and the long-sweep pass criterion (direct
          and twisted) leaves 0 of 2,020,802 windows for `r = 30..100`,
          `a = 3..8`.  FM-SEC27 (luna_max_venus) also noted that the pass
          criterion proves the four earlier grid windows
          (`fm39/w_pass_four_repro.py`).
        - Full scan, every row with `r <= 100` and `a <= 150` (14,798 rows,
          130,875,374 consumer windows): residual 0
          (`fm39/w_residual_rows3.py`).  So the three-`h` G0 branch is
          proved, for all labels, at every `r <= 100`, `a <= 150`.
          At every `r` the band G0B, the centered-window (equal-label)
          sector and the reduction G0E to (E) hold.
        - FM-SEC32 (luna_max_venus): no uniform exhaustiveness argument.
          The pointwise rotation bound `theta_k <= arccos rho_k` fails,
          e.g. at `(a,e,k) = (56,5,36)` (88.3 degrees against 34.7).  So
          the geometric criteria resist conversion into analytic ones.
      - Open: a uniform argument that the criteria are exhaustive.  The
        clearest target is the decreasing-side strip `a >= e + 2`,
        `x >= N/2`, `1 < u < C`.
    - *Recurrence model (main agent, `fm39/const_coeff_check.py`).*
      - The base rows satisfy
        `(k+1) c_(k+1) = (a-e) c_k - (N-k+1) c_(k-1)`.  So
        `g_(k+1) = M_k g_k`, with `M_k = [[alpha_k, -beta_k],[1,0]]`,
        `alpha_k = (a-e)/(k+1)` and `beta_k = (N-k+1)/(k+1)`.
      - Then `D_k = H_k(g_k)`, where
        `H_k(p,q) = p^2 - alpha_k p q + beta_k q^2`.  It is positive
        definite in the oscillatory zone, and `H_k(M_k v) = beta_k H_k(v)`.
      - If the coefficients were constant (`M` elliptic, rotation `omega`),
        then `D_(x+j) = beta^j H` and
        `T = beta^((n-1)/2) U_(n-1)(cos omega) H`, with `n = C+1`.  So W
        is `sum_(j<n) beta^j >= beta^((n-1)/2) |U_(n-1)|`, which holds by
        `|U_(n-1)| <= n` and AM--GM.
      - In the same model, (E) with `s = i-j` and `beta < 1` becomes
        `(1-beta^s) H >= |1-beta| beta^((s-1)/2) |U_(s-1)| H`, again true.
      - Checked in 560,000 windows.  So `D_k` is a discrete Sonin
        function, and Theorem OL is its monotonicity.  W and (E) are both
        comparisons of this kind with slowly varying coefficients, which
        points to a Sturm/Sonin comparison proof.
    - **W splits into two inequalities (main agent, `fm39/sonin_route.py`,
      `fm39/sonin_route2.py`).**  Put `T_C(x) = T(x, x+C)`, so `T_0 = D`,
      and `n = C+1`.
      - (i) `|T_C(x)| <= n sqrt(D_x D_(x+C))`.
      - (ii) `sum_(k=x..x+C) D_k >= n sqrt(D_x D_(x+C))`.
      - Together they give W.  On base rows (`r <= 8`, `a < 40`) both held
        in all 122,920 windows, of every length and every sweep.
      - Random rows:
        - (i) held on every real-rooted family: 22,251 windows with mixed
          roots, 24,124 with positive coefficients, 26,564 reciprocal
          (`|t| >= 2`).
        - (ii) held on the positive and reciprocal families, and failed
          once on a non-reciprocal mixed row (ratio 0.988, `C = 2`).
      - At `C = 1`, (i) is the real-rooted cubic discriminant inequality
        `4 D_k D_(k+1) >= T(k,k+1)^2`.  In the constant model it is
        `|U_C| <= C+1`.
      - Plucker recurrence (symbolic):
        `T_(C-1)(x) T_(C-1)(x+1) = D_x D_(x+C) + T_C(x) T_(C-2)(x+1)`, with
        `T_(-1) = 0`.  This is the analogue of
        `U_(C-1)^2 - U_C U_(C-2) = 1`.
      - (ii) follows from log-concavity of `D` along the window (AM--GM).
        But `D` is not log-concave on 93 of 280 base rows.  The failures
        are mild and spread over the flanks, e.g.
        `(1299, 2937, 6655)` at `r = 3`, `a = 11`.
      - For (E): `W = T_(s-1)(j) - T_(s-1)(j+1)` and
        `D_j - D_i = T_0(j) - T_0(i)`, with `s = i-j`.
    - **Theorem LD (inequality (i) for every real-rooted row; FM-SEC15,
      luna_max_venus; verified by the main agent, `fm39/ld_repro.py`,
      `fm39/ld_check.py`).**  For every real polynomial with only real
      roots and every support window,
      `|T_C(x)| <= (C+1) sqrt(D_x D_(x+C))`.
      - *Local Jensen polynomial.*  With `n = C+2` and `m = x+C+1`,
        `Q(t) = sum_(j=0..n) C(n,j) c_(x-1+j) t^j = [z^m] P(z)(z+t)^n` is
        real-rooted.  For `t_i` in the upper half-plane,
        `F = P(z) prod (z + t_i)` has all its `z`-roots in the closed lower
        half-plane.  By Gauss--Lucas so does `d^m F/dz^m`.  Hurwitz at
        `z -> 0` makes `d^m F/dz^m (0, t)` real stable or identically
        zero, and its diagonal `t_i = t` is real-rooted.  The
        identically-zero branch, e.g. `P = z^5`, `x = 0`, `C = 1`, makes
        the inequality trivial.  (Repair from FM-CHK36.)
      - FM-CHK36 (luna_max_jupiter): Theorem G0B ACCEPT on all five items,
        and it covers every support window, not only the consumer range.
        Theorem LD ACCEPT after the zero-branch repair above.
      - *Coefficient inequality.*  If `q = sum C(n,j) a_j t^j` is
        real-rooted, then
        `|a_1 a_(n-1) - a_0 a_n| <= (n-1) sqrt((a_1^2 - a_0 a_2)(a_(n-1)^2 - a_(n-2) a_n))`.
        Proof: by the root-sum identities, both factors and the left side
        are sums over pairs of roots, and Cauchy--Schwarz compares them.
        Degenerate cases follow by limits.
      - With `a_j = c_(x-1+j)` this is exactly (i).  `C = 1` is the cubic
        discriminant.
      - *Checks.*  The agent script reruns: 778,800 base windows, with the
        root-sum identities symbolic in degrees 2..8.  The main agent
        found every one of 3,870 random real-rooted windows gives a
        real-rooted `Q` (exact `real_roots`), and no violation of (i).
    - **What remains of W** (main agent, `fm39/w_residual.py`).
      - With the G0B energy bound `S >= max(D_x,D_y) + C min(D_x,D_y)`,
        Theorem LD gives `S >= |T|` whenever `u = sqrt(D_max/D_min)` lies
        outside `(1, C)`, since `u^2 + C >= (C+1) u` iff
        `(u-1)(u-C) >= 0`.
      - Consumer windows (`e` odd, `2x >= N-C`, `C >= 3`, `r <= 10`,
        `a < 60`), 233,865 in all:
        - `T <= 0`: 53.36%;
        - LD plus energy: 36.33%;
        - the G0B window region: 9.16%;
        - Theorem WS: 0.82%;
        - residual: 759 windows (0.32%).
      - Every residual window has large slack in (ii): the ratio
        `S/((C+1) sqrt(D_x D_y))` lies between 3.2 and 29.  So W, and with
        it the G0 branch, now reduces to (ii) on windows with
        `1 < u < C`.
      - Adding FM-SEC20's criteria (`e = 1`, and the reflection plateau)
        and the G0E rows (`min(a,e) <= 2`, strips) leaves 4 of the
        233,865 consumer windows (`fm39/w_residual2.py`,
        `fm39/w_residual3.py`).  All four have `r = 4`, `a = 56..59`,
        `C = 10, 11` and `u` just below `C`, and (ii)-ratios 3.3 to 4.0,
        so rows with `a` far above `e`.  Open: that the union of criteria
        is exhaustive beyond the grid.
    - **(E) at `q = 1` for exponent gaps 2 and 3 (FM-SEC17, luna_max_mars;
      verified by the main agent, `fm39/e_q1_repro.py`,
      `fm39/e_q1_route_kill.py`).**
      - For `a, e >= 3`, `|a-e| in {2,3}` and `i = j+2`, (E) holds.  This
        includes root-crossing pairs, e.g. `(a,e,j,i) = (5,3,5,7)` with
        `(L, W, slack) = (30, 10, 20)`.
      - *Proof.*  Fold to `a = e + d`.  The coefficient rows are explicit
        in `b_t = C(e,t)`, and `(L -+ W)/b_t^2` factors into products of
        nonnegative terms (four symbolic cases, `d = 2, 3`, `k` even or
        odd).
      - The outer endpoint `(j,i) = (N-1, N+1)` holds for every `a, e`
        with `N >= 2`.  There `D_(N-1) - D_(N+1) = (delta^2 + N)/2` and
        `|W| = |delta|`.
      - Kill (a proof route): allowing unrestricted recurrence starting
        values makes the `L - W` form indefinite (`N = 18`, `a-e = -12`,
        `j = 17`).  On the actual row, `L - W = 93`.
      - It covers fixed gaps only; `|a-e| >= 4` at `q = 1`, and `q >= 2`
        across roots, remain open.
      - *Row class of (E) (main agent, `fm39/e_class.py`).*  On
        reciprocal rows `(1+z)^a (1-z)^e` times one extra factor, (E)
        held:
        - with real reciprocal pairs (`|t| >= 2`): 5,714 pairs, 0
          violations;
        - with complex root quadruples whose circle weight
          `((x+al)^2 + be^2)` is log-concave on `[-2,2]`: 6,741 pairs, 0
          violations.
        It failed:
        - without log-concavity: 637 violations in 248 of 400 rows;
        - with unit-circle roots: 1,377 violations.
        So (E) lives on the class of reciprocal rows with log-concave
        circle weight, which is wider than real-rooted.  The three-factor
        value tracks the same class.  This points to an analytic
        (Prekopa--Leindler type) mechanism.
      - *Fourier representation (main agent, `fm39/fourier_rep.py`,
        machine precision, both parities).*  Let `F` be the real circle
        weight of the row (a cosine series for even `e`, a sine series
        for odd `e`), `m_k = k - N/2`, `u = t - p`, `v = t + p`.  Then
        - `D_k = <F(t)F(p) cos(m_k u)(1 - cos v)>`;
        - `T(x,x+C) = <F(t)F(p) cos(M u) 2 sin((C+1)v/2) sin(v/2)>`, with
          `M = x + C/2 - N/2`.
        For base rows `F(t)F(p)` is proportional to
        `K(2 cos(u/2), 2 cos(v/2))`, the kernel itself, so this is a
        Chebyshev-T re-expression.  In it, (E) compares `q+1`
        consecutive unit decrements in `M` with `|U_q(cos(v/2))|` times
        the middle one, averaged in `v`.  The pointwise version already
        fails at `a = e = 1`, so a proof must use the `v`-average.
      - FM-SEC26 (luna_max_uranus): Chebyshev-difference form of (E).
        With `q = i-j-1` and `eta = (i+j-N)/2`,
        `D_j - D_i -+ W_ij = <4K sin(eta u) sin(u/2) sin^2(v/2) [U_q(cos(u/2)) -+ U_q(cos(v/2))]>`.
        The weight changes sign for `q >= 1`, so log-concavity alone does
        not give a positive comparison.  No movement.
      - FM-SEC22 (luna_max_mars): no uniform proof.
        - Exact drift identity:
          `S(j,j+2) = d D_j/(j+1) + R_j`, with
          `R_j = c_(j-1)(d c_(j+1) - (N+2) c_j)/((j+1)(j+2))`.  So (E) at
          `q = 1` reads
          `D_j - D_(j+2) >= |d (D_j/(j+1) - D_(j+1)/(j+2)) + R_j - R_(j+1)|`.
          Here `R_j - R_(j+1)` is the recurrence drift that the
          constant-coefficient model omits.
        - Screen: 33,200 `q = 1` pairs (`a, e <= 40`), no failure.  The
          least open-slice slack is 4.
        - Demand sheet re-derived at this stall: the consumer uses (E) with
          both signs at both parities of `e` (`h h` and `hat S hat S` at
          even `e`; `h hat S` at odd `e`).  So the absolute-value form is
          exactly what is used.
    - **(E) on an energy-drop region, uniformly (FM-SEC25, luna_max_mars;
      verified by the main agent, `fm39/e_ld_region_repro.py`).**
      - Since `W_ij = S(j,i) - S(j+1,i+1)`, and both chords are long Turan
        determinants of window length `s-1` (`s = i-j`), Theorem LD and
        Theorem OL give (E) whenever
        `D_j - D_i >= s (sqrt(D_j D_(i-1)) + sqrt(D_(j+1) D_i))`.
        A simpler sufficient condition is
        `D_j/D_(i-1) >= (s + sqrt(s^2+1))^2`.
      - This is uniform in `a`, `e` and `q`: 298,234 of 438,221 pairs
        (`a, e <= 40`) lie in the region, with no failure.
      - Outside it lies e.g. the root-crossing witness `(3,5,5,8)`.
    - **Residual of (E)** (main agent, `fm39/e_residual.py`; RF omitted,
      so this is conservative).  Over 438,221 pairs:
      - LD energy drop: 58.27%;
      - strips or `min(a,e) <= 2`: 12.90%;
      - `q = 0`: 6.77%;
      - `q = 1`, gaps 2 and 3: 0.69%;
      - outer endpoint: 0.27%;
      - residual: 92,468 pairs (21.10%), spread over every `q` and every
        gap.
    - **(E) as a polygon inequality for a second curve (main agent,
      `fm39/psi_curve.py`, `fm39/psi_split.py`).**
      - Put `psi_k = (c_k, B_k)`, `B_k = c_(k-1) + c_(k+1)`.  Then
        `W_ij = psi_j ^ psi_i` and
        `delta_k := D_k - D_(k+1) = psi_k ^ psi_(k+1)`, which is `>= 0` by
        Theorem OL.  So (E) is exactly
        `|psi_j ^ psi_i| <= sum_(k=j..i-1) psi_k ^ psi_(k+1)`, the
        chord-versus-area inequality of Theorem WS and W, for the curve
        `psi`.
      - Left turns of `psi` are (E) at `q = 1` with the minus sign:
        `(psi_(k+1)-psi_k) ^ (psi_(k+2)-psi_(k+1)) = D_k - D_(k+2) - W_(k+2,k)`.
        So the convex-polygon argument gives (E) at every `q`, both signs,
        on `psi`-windows sweeping `< pi`, from (E) at `q = 1`.  But only
        8% of pairs sweep `< pi`; 74% sweep `> 2 pi`.
      - *Sibling rows (FM-SEC28, luna_max_mars).*  With
        `F = (1+z)^(a+2)(1-z)^e` and `G = (1+z)^a (1-z)^(e+2)`, one has
        `F - G = 4zP` and `F + G = 2(1+z^2)P`.  Hence
        `W_ij = (f_(j+1) g_(i+1) - g_(j+1) f_(i+1))/4`, and `delta_k` is the
        adjacent mixed determinant.  So (i)-psi is an LD-type bound for
        the mixed curve of two real-rooted rows sharing the factor `1+z`
        (for `a >= 1`).  The `a = 0` failures show that real-rootedness of
        the two rows alone is not enough.  Reproducer:
        `fm39/psi_sibling_repro.py`.
      - Split, parallel to W's (i)+(ii):
        - (i)-psi: `|W_ij| <= (i-j) sqrt(delta_j delta_(i-1))`.  It held
          in 402,926 of 403,340 pairs with `q >= 1`; the failures are at
          `a = 0`, which is already proved.
        - (ii)-psi: `sum_(k=j..i-1) delta_k >= (i-j) sqrt(delta_j delta_(i-1))`.
          It held in 400,915.
        - Both held in 349,832 of 352,038 open-region pairs (99.4%).
          (ii)-psi fails e.g. at `(4,6,5,8)`.
        - Correction (FM-SEC29, luna_max_jupiter): (i)-psi also fails at
          a few even-`e` pairs, e.g. `(4,36,20,22)` and its folds, not only
          at `a = 0`.  Odd-`e` consumer pairs are unaffected.  FM-SEC29
          proves (ii)-psi on several uniform subregions (span 2,
          endpoint-imbalance, locally log-concave increments, `e = 1`).
          Its odd-`e` open residual is 991 grid pairs, all positive with
          least (E) slack 142.  Reproducer: `fm39/psi_ii_repro.py`.
    - **Metric chords shrink the (E) residual (main agent,
      `fm39/e_metric.py`).**  Bound each chord `T(j,i-1)` and `T(j+1,i)`
      of `W_ij` by the smaller of Theorem LD and the FM-MECH25 metric
      bound `sqrt(K D_x D_y / ((2 sqrt V - R_+)(2 sqrt V - R_-)))`.  Then
      (E) holds whenever `D_j - D_i` exceeds the sum.  This proves
      another 48,872 pairs.
      - Residual of (E) on `a, e <= 40`: 43,596 of 438,221 pairs (9.95%,
        down from 21.1%).
      - The residual concentrates at the shortest gaps: `q = 1` 17,944,
        `q = 2` 11,970, `q = 3` 6,908, `q = 4` 3,948, `q = 5` 1,988,
        `q >= 6` 838.
      - At `q = 1` the chords are single steps and separate absolute
        bounds lose the cancellation.  So `q = 1` across all gaps is the
        sharpest open core of (E).
    - **Theorem EQ1: (E) at `q = 1` for all `a, e`, from Theorem OL on the
      neighbouring rows (main agent, `fm39/e_q1_ol.py`).**
      - Put `psi_k = (c_k, B_k)`.  Then `psi_(k+1) - psi_k` is the
        `psi`-vector of `f = (1-z)P` at `k+1`, and `psi_k + psi_(k+1)` is
        that of `s = (1+z)P`.  Hence, exactly:
        - `D_j - D_(j+2) - W_(j+2,j) = D^f_(j+1) - D^f_(j+2)`, with
          `f = (1+z)^a (1-z)^(e+1)`;
        - `D_j - D_(j+2) + W_(j+2,j) = D^s_(j+1) - D^s_(j+2)`, with
          `s = (1+z)^(a+1)(1-z)^e`.
      - Both rows are base rows of degree `N+1`, and `j >= N/2` gives
        `2(j+1) > N+1`.  So Theorem OL makes both sides `>= 0`.  QED.
        This replaces the fixed-gap factorizations of FM-SEC17 and proves
        `q = 1` for every gap.  Checked on all 33,200 `q = 1` pairs
        (`a, e <= 40`).
    - **The multiplier identity for (E) (main agent,
      `fm39/e_multiplier_symbolic.py`, `fm39/e_ol_lp.py`,
      `fm39/e_via_olcyc.py`).**  With `delta_k[R] = D_k(RP) - D_(k+1)(RP)`,
      - `D_j - D_i - W_ij = sum_(i0 < q/2) delta_(j+q)[z^i0 - z^(q-i0)]`;
      - `D_j - D_i + W_ij = sum_(i0 < q/2) delta_(j+q)[z^i0 + z^(q-i0)]`,
        plus `2 delta_(j+q)[z^(q/2)]` when `q` is even.
      This is symbolic for `q = 1..7`.  For general `q` (FM-CHK39,
      luna_max_jupiter), put `u = j+r`, `v = j+q-r`, `m = q-2r` for
      `r < q/2`.  Then
      `delta_(j+q)[z^r (1 - sigma z^m) P] = (psi_v - sigma psi_u) ^ (psi_(v+1) - sigma psi_(u+1))`.
      Summing, the adjacent areas give every `delta_j..delta_(j+q)`
      (except the middle one when `q` is even), and the cross terms
      cancel in pairs to `W_ij`.  Each term lies beyond its row's OL
      centre, since `2(j+q-r) - (N+m) = 2j + q - N >= q > 0`.
      - Each term is a one-label (OL) difference of the row
        `(1 -+ z^m) P`, `m = q - 2 i0`, beyond its centre.  So (E) follows
        from "Theorem OL for cyclotomic multiples `(1 -+ z^m)P`", the
        m = 1, 2 cases of which are base rows.
      - Consequence: (E) at `q = 2` with the minus-W sign is Theorem OL for
        `(1 - z^2)P`, so it is proved for all `a, e`.
      - Termwise OL for the cyclotomic multiples holds in 686,420 of
        704,076 open-region (pair, sign) cases (97.5%).  In the rest a
        negative term is outweighed within the sum, e.g. `[-10, 40]`.
        Termwise OL for `m >= 3` is unproved: Theorem OL's proof is
        Krawtchouk-specific.
      - Kills of linear routes (knob: linear certificates, which the
        consumer does not require).
        - The `q = 2` plus sign is not a nonnegative combination of OL
          differences and proved `q = 1` or `q = 2` (minus) forms of
          base-family rows (`fm39/e_recursive_lp.py`, multipliers of
          degree `<= 6`).  Nor is `q = 3`.
        - In word terms the row `(1 -+ z^m)P` inserts the virtual factor
          `h_(m-1) - h_(m-3)` (or `hat S_m - hat S_(m-2)`).  The identity is
          the Clebsch--Gordan telescoping of `U_q`, so termwise positivity
          is not expected; the sum is what the consumer uses.
      - What remains at `q = 2` is the plus sign,
        `delta_(j+2)[(1+z^2)P] + 2 delta_(j+1)[P] >= 0`.
        - FM-SEC35 (luna_max_venus) restated it via `(1+z^2)P = (F+G)/2`,
          with `F = (1+z)^2 P` and `G = (1-z)^2 P`, as
          `delta_(j+2)[F] + delta_(j+2)[G] >= 4 delta_(j+1)[P]`.  All
          three rows are base rows; open for `a >= e+2`, `e >= 3`.
        - The Cauchy--Schwarz bound on the cross term fails at
          `(a,e,j) = (4,6,5)`.  No proof yet (FM-SEC38 is applying
          Theorem OL's proof method).
    - **Uniform metric regions for (E) (FM-MECH26, astra_max_ceres;
      reproducer `fm39/e_metric_mech26_repro.py`).**
      - Put `sigma = N+2`, `A_k = c_(k-1) - c_(k+1)` and
        `v_k = (c_k, B_k/2)`.  Then `d c_k - sigma B_k/2 = -X_k A_k/2` and
        `D_k = c_k^2 - B_k^2/4 + A_k^2/4`.
      - For `0 < t < 4V` the form `Q_t(f,g) = f^2 - g^2 + (d f - sigma g)^2/t`
        is positive definite, and `Q_t(v_k) = F_k(t)/(4t)` with
        `F_k(t) = 4t D_k + (X_k^2 - t) A_k^2`.  Since
        `W_ij = 2 det(v_j, v_i)`, one metric bounds the whole W,
        keeping the two-chord cancellation:
        `4t(4V-t) W_ij^2 <= F_j(t) F_i(t)`.
      - Taking `t = X_i^2` gives region `(M_i)`, and `t = X_j^2` region
        `(M_j)`.  Each proves (E), both signs, for all `a, e`.
      - It gives an independent proof of all `q = 1` pairs via an OL
        identity (agreeing with Theorem EQ1), and makes the strict
        short-`psi`-arc argument unconditional.
      - The union closes 92,466 of the 92,468 earlier residual grid pairs.
        Global coverage is not claimed.
    - **Binary-form method for (E) (main agent, `fm39/e_q2_quadratic.py`,
      `fm39/e_binary_form.py`, `fm39/e_residual3.py`).**
      - The three-term recurrence (valid with zero extension at both ends)
        writes every coefficient of the window `j-1..i+1` through
        `(c_j, c_(j+1))`.  So each sign of (E) at `(j, i)` is an explicit
        binary quadratic form in `(c_j, c_(j+1))`, with coefficients
        rational in `(d, N, j)` and positive denominators.  If it is
        definite (discriminant `< 0`, leading coefficient `> 0`), (E)
        holds with no information on the ratio.  This is step 1 of
        Theorem OL's proof, extended to every `q`.
      - Definite fraction in the open region (`a, e <= 40`), `q = 1..6`,
        both signs: from 96.5% (`q = 1`) up to 99.8% (`q = 6`), rising
        with `q`.
      - **Combined, (E) holds on every pair with `a, e <= 40`**:
        - LD energy drop or metric chords: 68.02%;
        - strips, `min(a,e) <= 2`, `q <= 1` (Theorem EQ1), outer
          endpoint: 26.13%;
        - definite binary forms: 5.85%;
        - residual: 0 of 438,221.
        The binary-form tests are exact rational arithmetic.  The LD and
        metric coverage tests use floating point with relative tolerance
        `1e-12`.  FM-CHK39 certified all 298,072 float-accepted decisions
        on this grid by exact integer comparisons, with rationalized
        square-root bounds.
      - FM-CHK39 (luna_max_jupiter), own exact code: ACCEPT for Theorem
        EQ1 (every boundary case), the multiplier identity (all-`q`
        proof above), and the binary-form method (the recurrence holds at
        `k = 0` and `k = N`, denominators are positive, and definiteness
        implies the value at the actual seed).  No non-definite
        three-factor form was found on rerun.
      - Systematic scan, every row with `3 <= a, e <= 80`: 6,329,895
        pairs, residual 0 (`fm39/e_residual_rows.py`).  Rows with
        `min(a,e) <= 2` or `|a-e| <= 1` are covered by earlier theorems.
      - Random scan, 600 rows with `a, e <= 300`: 7,885,359 pairs,
        residual 0, with Theorem RF used on 113 pairs
        (`fm39/e_residual_sample2.py`; D rescaled per row before the
        floating-point LD/metric tests).  An earlier 150-row sample without
        RF left 10 pairs in rows `(180,3)` and `(188,5)` (small `e`), all
        of which RF covers.
      - So the two-label stratum (all three sign patterns) holds at every
        level for all `a, e <= 80`, and through Theorem G0E so does the G0
        branch, whose rows have `e` odd.
      - Exactness: all 5,274,886 floating-point LD/metric decisions of the
        `a, e <= 80` scan are certified by exact integer comparisons
        (`fm39/e_cert_rows80.py`, FM-CHK39's certifier).
      - **Uniform regions alone** (`fm39/e_uniform_residual.py`,
        `fm39/e_uniform_sample.py`).  Use only criteria stated for all
        `a, e`: the strips, `min(a,e) <= 2`, `q <= 1` (EQ1), the outer
        endpoint, the LD energy drop, metric chords, FM-MECH26's `(M_i)`
        and `(M_j)`, and Theorem RF.
        - On `3 <= a, e <= 40` they leave 0 pairs.
        - On 600 random rows with `a, e <= 300` (7,176,773 pairs) they
          leave 70 pairs, every one with `q = 2`.
        - The 70 are all lopsided rows (`min(a,e)` in 3..11 (FM-SEC42),
          `max(a,e)` up to 293) with `W > 0`, where the plus sign is
          trivial: `D_j - D_i + W >= D_j - D_i >= 0` by OL.  At `q = 2`
          the minus sign is OL on `(1-z^2)P`, so every `q = 2` pair with
          `W >= 0` is proved.
        - Adding that criterion (`fm39/e_uniform_residual2.py`,
          `fm39/e_uniform_sample2.py`) gives residual 0 on
          `3 <= a, e <= 40` and on all 7,141,576 sampled pairs up to
          `a, e <= 300`.
        - So, on everything tested, (E) is covered by criteria each proved
          for all `a, e`.  What remains is exhaustiveness: that this union
          of explicit inequalities covers every pair.  The finite-`q`
          route is closed as posed (FM-SEC41 below); the `q = 2`, `W < 0`
          plus sign is FM-SEC38, FM-SEC44.
      - **FM-SEC41 (luna_max_jupiter): no finite `q0` for this menu.**
        - Family: `e = 3`, `X_k = 2k - N`; `x_n` is the largest lattice
          point below `sqrt(3N+4)` (the positive root of `K_3(X; N+2)`),
          `y_n` the smallest with `16 y_n^2 >= 75N`, and
          `(j, i) = ((N+x_n)/2, (N+y_n)/2)`.  Then
          `q_n = (sqrt 3/8) sqrt N + O(1)`.
        - Exact samples `(a, j, i, q)`: `(2000,1040,1050,9)`,
          `(4000,2056,2070,13)`, `(8000,4078,4099,20)`,
          `(16000,8111,8139,27)`.  LD, `(M_i)` and `(M_j)` fail (exact
          integers, rerun by the main agent).  The window straddles a
          Krawtchouk root, so RF does not apply.  (E) holds, with limiting
          normalized slack `> 0.1664` (rational exponential bounds).
        - Grid replay `a, e <= 60` (2,104,931 pairs): exact LD with
          `(M_i)`/`(M_j)` leaves 232 pairs, all `q = 1` (EQ1).
      - **Short `psi`-arcs cover the family (main agent,
        `fm39/psi_arc_family.py`).**  All four samples have
        `W_(j,k) > 0` for every `j < k <= i` (sweeps 3.086, 3.106, 3.134,
        3.055 rad).  The polar angle of `psi` is monotone (Theorem OL), so
        this is sweep `< pi`, and the convex-polygon argument from EQ1
        gives (E), both signs.
        - Membership is decided by the sign of `W_(j,k)` at the actual
          row, as for `q = 2, W >= 0`.  So it is an explicit criterion,
          not a region in `(a, e, j, i)`.
      - **Extended menu (main agent, `fm39/e_extreme_residual.py`, exact
        integers).**  Menu: `q <= 1`, outer endpoint, short `psi`-arc,
        `q = 2` with `W >= 0`, LD, `(M_i)`, `(M_j)`, FM-MECH26's general
        bound `4t(4V-t) W^2 <= F_j(t) F_i(t)` at `t = 2V`, RF.  No binary
        forms and no float chords.
        - Lopsided rows `a in {500, 1000, 2000}`,
          `e in {3,...,8,10,12,16,20,30}`: 7,351,913 pairs, residual 0
          (short arcs 78%, LD 21%).
      - **Fixed-`e` limit of (E) and of the menu (main agent,
        `fm39/e_continuum_limit.py`, `fm39/e_continuum_cover.py`,
        `fm39/e_center_mt.py`).**
        - With `X_k = t sqrt N` and `a -> infinity`, the row tends to
          `F = D^e exp(-t^2/2)`, which solves `F'' + t F' + (e+1) F = 0`.
          Then `D_k/S -> G = F'^2 + t F F' + (e+1) F^2` and
          `W_ij/S -> th F'(th) F(et) - et F(th) F'(et)`.  Checked against
          exact rows at `a = 4000, 16000` (3 decimals).
        - So (E) tends to
          `(E_inf,e)`: `G(th) - G(et) >= |th F'(th) F(et) - et F(th) F'(et)|`
          for `0 <= th < et`.  This is the chord-versus-area inequality
          for the curve `(F, -t F')`, the sheared limit of `psi`.
        - Numerically `(E_inf,e)` holds for `e <= 40` and `e = 60, 80,
          120`.  The ratio `|W| / (G(th) - G(et))` tends to 1 only on the
          diagonal, where the minus sign has zero first-order slack (the
          limit of EQ1).
        - Limit of the menu.  LD drops out, since its factor `i - j` grows
          like `sqrt N`.  Short arcs, `(M_i)`, `(M_j)` and RF cover every
          grid point of `0 <= th < et <= 2 sqrt(e+1) + 3` (`e <= 30`)
          except `th = 0`, `et >= 2 sqrt(e+1)`, `e` even.  There `(M_i)`
          needs `X_i^2 < 4V`, and `(M_j)` needs `X_j != 0`.  (For odd `e`,
          `c_(N/2) = B_(N/2) = 0`, so `W = 0` at `j = N/2`.)
        - The general bound at `t = 2V` closes that line: in the limit for
          every even `e <= 40`, and exactly at `j = N/2` on rows
          `(500,4)`, `(1000,4)`, `(2000,6)`, `(2000,10)`, `(4000,4)`,
          `(4000,20)` (5,964 pairs, no failure).
        - So for fixed `e` no family escapes the extended menu in the
          limit.  This is evidence, not a proof.  A proof needs error
          bounds uniform in `a`, and the proportional regime `e ~ N` is a
          different limit.
      - **Two-case split of (E) and Conjecture HT (main agent,
        `fm39/e_sweep_split.py`, `fm39/e_sweep_split_t.py`, exact
        integers).**
        - Rows with `a >= e+2`, `e >= 3` (the rest is proved or follows by
          `a <-> e`), pairs with `q >= 2`.  Case (b): the `psi`-arc from
          `j` to `i` sweeps `< pi`; proved by the convex-polygon argument
          from EQ1.  Case (c): it sweeps `>= pi`.
        - On every case-(c) pair tested, FM-MECH26's general bound closes
          (E):
          `(HT)`: some `t in (0, 4V)` has
          `4t(4V-t)(D_j - D_i)^2 >= F_j(t) F_i(t)`.
          - Grid `3 <= e <= a-2 <= 78`: 2,638,251 pairs, no exception.
          - Grid `3 <= e <= a-2 <= 198` (19,306 rows): 108,537,108 pairs,
            no exception, each closed at `t in {X_i^2, X_j^2, 2V}`.
          - 24 large rows, up to `(2000,20)`, `(500,498)`, `(1000,500)`:
            1,083,099 pairs, no exception.
          - With `a <= 60`, `t = 2V` alone covers 803,778 of the 812,334
            pairs, and `t = X_i^2` or `t = X_j^2` covers the rest.  `t = 2V`
            is needed only at `X_j = 0`.
        - Since `4t(4V-t) W^2 <= F_j(t) F_i(t)` is proved, (HT) on case
          (c), with EQ1 and the convex-polygon argument, proves (E) for all
          `a, e`.
        - FM-CHK41 (luna_max_saturn, own exact code).
          - ACCEPT: the single-metric bound for every `t in (0,4V)`,
            including `k = N+1` (zero extension: `(c,B,A,D) = (0,1,1,0)`)
            and `X_k = 0`.
          - REPAIR, open: the short-arc convexity argument needs strict
            steps.  Weak conditions allow flat steps that backtrack along a
            ray, e.g. the point chain `(1,0),(3,3),(1,1),(2,2),(4,4)`, whose
            doubled area is `-1`.  This is not a row counterexample.  What
            is needed is strict OL (`delta_k > 0` for `2k > N` whenever
            `psi_k != 0`) or a flat-block lemma.  An isolated flat step with
            strict neighbours is forced by EQ1 to be a repeated point, which
            is harmless.  Assigned as FM-SEC64.
          - REPAIR, applied.  Zero anchors form their own case: `psi_k = 0`
            only at `N = 2k` with `e` odd, where `W_(k,i) = 0` and (E) is
            OL.  The centre sign at even `N = 2m`:
            `delta_m = c_m^2 (4(m+1)^2 - (a-e)^2)/(2(m+1)^2(m+2)) >= 0`
            for even `e`, and `delta_m = 0` for odd `e`.  The fixed-`e`
            normalization is `C_N = C(N, floor(N/2)) N^(-e/2)`,
            `S = 4 C_N^2 / N`.
          - Conclusion: the split (strips, fold, `q <= 1`, short arcs,
            case (c)) is complete, conditional on (HT) and on the
            strictness repair.
          - FM-SEC64 (luna_max_saturn): the strictness repair is done, by a
            flat-block lemma (`fm39/sec64_flatblock_repro.py`, resultant
            rerun exactly).
            - A zero `c_k` never gives a flat step in the upper half:
              `delta_k = p(N-k+1) c_(k-1)^2/(k+1)^2 > 0`.
            - There are never two consecutive flat steps.  The resultant
              of the two flatness conditions is
              `(p+2)^2 (k+2)^2 (N+2-d)(N+4-d)(N+2+d)(N+4+d) > 0`.
            - By EQ1, an interior flat step is a repeated vertex and a
              leading flat edge moves outward.  A trailing inward flat edge
              (main agent) is also harmless, since
              `W_ij = lambda W_(j,i-1)` with `lambda <= 1`.
            - So the short-arc case of (E) is complete: (E) now rests
              only on (HT) for case (c).
            - Strict OL itself is proved for `e <= 3` and `3p >= N-2`, and
              screened to `a, e <= 60`.  The band boundary `F = 0` has
              infinitely many Pell-type tuples, e.g. `(e,N,p) = (2,49,1)`;
              those in the box `e <= 100`, `N <= 1500` have ratio `!=` the
              double root.  The consumer does not need that audit.  That gives the two-label stratum and, through Theorem
          G0E, the whole G0 branch at every level.  LD, RF and the binary
          forms are then not needed.
        - Heuristic (unproved).  `det Q_t = (4V-t)/t`.  If `D` behaved
          like the square norm in one fixed metric, then over a half-turn
          of `psi` it would shrink by `exp(-pi sqrt(kappa))`, with
          `kappa = 4t/(4V-t)`.  At `t = X_i^2`, (HT) needs only
          `D_i/D_j <= rho*`, where `(1-rho*)^2 = kappa rho*`: about
          `1 - sqrt(kappa)` for small `kappa` and `1/kappa` for large
          `kappa`.  That leaves a factor `pi` of room in the exponent.
        - FM-SEC47 (luna_max_jupiter, continuum route) and FM-SEC48
          (luna_max_neptune, discrete route) are attacking (HT).
        - Margin map.  The true margin is `(2V(D_j - D_i)/perm)^2`, with
          `perm = |m_j n_i| + |n_j m_i|` (closed form below;
          `fm39/e_ht_perm_margin.py`).
          - Over all 1,528,422 case-(c) pairs with `a <= 70` the minimum
            is 2.783, at `(70,68,69,72)` (`q = 2`, `X_j = 0`, even `e`).
            Off the centre line it is 4.63 for even `e` and 4.02 for odd
            `e`.
          - Along `a = e+2` the centre-line minimum tends to `25/9`
            (2.7779 at `(400,398)`).  On large lopsided rows up to
            `(1000,3)` and `(500,250)` it is at least 4.25.
          - Correction (FM-SEC54, luna_max_jupiter).  An earlier value,
            `1.674` from `fm39/e_ht_margin.py`, searched only finitely
            many `t` and missed the limit `t -> 0` on the centre line.  It
            understated the margin; the true value at `(14,4,9,12)` is
            `(2843/132)^2`.
          So (HT) holds with room.  The room is smallest on the centre
          line near the diagonal, at `q = 2`.
        - `t = 2V` fails only on bulk pairs near the centre (all 4,271
          failures with `a <= 50` have `X_j^2 <= X_i^2 <= 2V`); there
          `t = X_i^2` or `X_j^2` works (`fm39/e_ht_2v_fail.py`).
          `Q_(2V)` has determinant 1, so `(HT)` at `t = 2V` reads
          `D_j - D_i >= 2 sqrt(Dt_j Dt_i)`, with
          `Dt_k = Q_(2V)(v_k) = D_k + (X_k^2 - 2V) A_k^2/(8V)`.
        - FM-SEC48 (luna_max_neptune): exact reduction.  With
          `P_k = 4D_k - A_k^2`, `Q_k = X_k^2 A_k^2`, `E = D_j - D_i`,
          `(HT)` holds iff `L = 4E^2 + P_j P_i > 0`, `0 < b < 8LV` and
          `b^2 >= 4 L Q_j Q_i`, where `b = 16VE^2 - P_j Q_i - P_i Q_j`
          (or `L = b = Q_j Q_i = 0`).  At `X_j = 0` the last condition is
          void.  Boundary families pass exactly.
        - Killed routes (knob: auxiliary monotonicity that `(HT)` does
          not use).
          - The fixed-metric radial drop.  `R_k = Q_t(v_k)` is not
            monotone along case-(c) arcs, and `R_i > R_j` occurs
            (`fm39/e_spiral_probe.py`).  One step with `delta_k = 0` has
            `R_k - R_(k+1) = -421/12` at `(a,e,k,t) = (5,3,4,48)`
            (FM-SEC48).  So the log-spiral heuristic above is withdrawn.
          - Monotonicity of the margin in `i`: `rho(j,i)` decreases in
            19,304 of 152,128 steps with `a <= 40`, e.g. `(6,3,5,9)`
            (`fm39/e_ht_monotone.py`).
      - **Closed form of (HT): the permanent inequality (main agent,
        `fm39/e_ht_mn.py`, exact on every row `a, e <= 30`).**
        - Put `m_k = sigma c_k - d B_k/2`, `n_k = X_k A_k/2` (the
          recurrence gives `d c_k - sigma B_k/2 = -X_k A_k/2`), and
          `p_k = (m_k, n_k)`.  Then:
          - `4V D_k = m_k^2 + mu_k n_k^2`, with `mu_k = 4V/X_k^2 - 1`.  So
            `D` is a diagonal metric in `(m, n)` whose weight decreases in
            `k` for `X_k > 0` and changes sign at the turning point
            `X_k^2 = 4V`.
          - `2V (D_k - D_(k+1)) = det(p_k, p_(k+1))`: the fan area of the
            `p`-curve equals the drop of its own varying norm.
          - `2V W_ij = m_j n_i - n_j m_i`, so (E) is
            `2V (D_j - D_i) >= |m_j n_i - n_j m_i|`.
        - `Q_t = (sigma f - d g)^2/(4V) + (d f - sigma g)^2 (1/t - 1/(4V))`,
          so the family is all positive combinations of two fixed squares.
          Minimizing the metric bound over `t` in closed form gives
          `(HT*)`: `2V (D_j - D_i) >= |m_j n_i| + |n_j m_i|`.
          This is (E) with the permanent of absolute values in place of
          the determinant.  It is (HT) together with its limits
          `t -> 0, 4V`.  On `a, e <= 30` the two differ only in 404
          equality cases with one product zero (for example `i = j+1`).
        - `(HT*)` is (E) together with
          `(E^R)`: `2V (D_j - D_i) >= |m_j n_i + n_j m_i|`.
        - Constant-metric model.  With `mu` fixed, the drop identity makes
          `p` a log-spiral in the metric `diag(1, mu)`, and over a
          half-turn `(HT*)` holds with a factor `pi` of room
          (`sinh(pi/sqrt(mu)) >= 1/sqrt(mu)`).  The right radius is
          `4V D_k`, which is monotone by OL.  The fixed-`Q_t` radius
          killed above is not.  FM-SEC54 is making this exact with the
          varying weight `mu_k`.
      - **Drop lemma and a rigorous partial chain (main agent,
        `fm39/e_ht_drop.py`, `fm39/e_ht_sinsum.py`).**  Take bulk
        case-(c) pairs `0 < X_j < X_i`, `X_i^2 < 4V`, so
        `mu_j >= mu_k >= mu_i > 0`, and put `rho_k^2 = 4V D_k`.
        - `|m_j n_i| + |n_j m_i| <= rho_j rho_i / sqrt(mu_i)`.  So `(HT*)`
          follows from the drop lemma
          `(DL)`: `D_i/D_j <= r*(mu_i)^2`, where
          `r*(mu) = sqrt(1 + 1/mu) - 1/sqrt(mu)`.
          `(DL)` involves only `D` at the two ends and the half-turn.  It
          holds on all 124,328 bulk case-(c) pairs with `a <= 40`.
        - Per step, exactly: in the metric `diag(1, M)` with
          `M >= mu_k, mu_(k+1)`, `2V delta_k = det(p_k, p_(k+1))` gives
          `sinh(x_k - x_(k+1)) >= sin(dth_k^(M)) / sqrt(M)`, with
          `x = log rho` and `dth^(M)` the elliptic angle step.  `sinh` is
          superadditive and `x` decreases (OL).  So `(HT*)` holds whenever
          `(b)`: `sum_k sin(dth_k^(mu_k)) / sqrt(mu_k) >= 1/sqrt(mu_i)`
          (local metric at each step).  This chain is a proof for each
          pair that meets `(b)`.
        - `(b)` holds on 1,297,485 of the bulk case-(c) pairs with
          `a <= 70`.  It fails on a thin set that persists at large `N`,
          e.g. `(400,3,202..206,241)` (ratio 0.75).  The failures have
          `X_i` near the turning point, where `mu_i -> 0` and the bound
          `|n_i| <= rho_i / sqrt(mu_i)` degenerates.  There `t = 2V`
          already closes `(HT)` on every tested pair.  So a proof can
          split: `(b)`-type estimates in the bulk away from the turning
          point, and the `t = 2V` bound near and past it.
        - The split works on everything tested
          (`fm39/e_ht_regions.py`, exact rows, float angles).  It covers
          all 800,154 case-(c) pairs with `a <= 60` and 8 large rows up to
          `(800,3)`, `(800,400)`:
          - inner, `X_i^2 <= 2V`: `(b)` holds, or `(b')` holds: the same
            chain with the exact ratio `perm/(rho_j rho_i)` in place of
            `1/sqrt(mu_i)` (12 pairs, all `q = 2`).  `(b')` is also a proof
            for its pair.  The centre `X_j = 0` has `n_j = 0`, and the first
            step uses the metric `mu_(j+1)`.
          - outer, `X_i^2 > 2V`: `(HT)` at `t = 2V`, i.e.
            `D_j - D_i >= 2 sqrt(Dt_j Dt_i)`.
        - At `a <= 120` (6,786 rows): all 13,634,429 case-(c) pairs pass,
          7,797,535 outer by `t = 2V` and 5,836,894 inner by `(b)`, 34 of
          them by `(b')`.
        - So (E) reduces to two statements on case-(c) pairs.
          - `(I)` (inner): `sum_(k=j..i-1) sin(dth_k^(mu_k)) / sqrt(mu_k)
            >= (|m_j n_i| + |n_j m_i|) / (rho_j rho_i)`.  This is a discrete
            phase-integral bound over a half-turn, and each term is an
            explicit function of the consecutive ratio `c_(k+1)/c_k`.
          - `(O)` (outer): `D_j - D_i >= 2 sqrt(Dt_j Dt_i)`.  Since
            `Q_(2V)(v) = ((sigma f - d g)^2 + (d f - sigma g)^2)/(4V)`,
            one has `4V Dt_k = m_k^2 + n_k^2 = |p_k|^2` (checked at
            30,752 points).  So `(O)` reads
            `2V (D_j - D_i) = sum_k det(p_k, p_(k+1)) >= |p_j| |p_i|`:
            the Euclidean fan area of the `p`-curve over an outer
            half-turn is at least the product of the end radii.
            FM-SEC55 (luna_max_venus) is attacking it.
            FM-SEC55 result: conditional.  Exact integer screens pass on
            777,100 outer pairs: 471,319 with `a <= 60` (least ratio
            `123904/23953`), large rows to `(2000,3)`, and `a = e+2` to
            `e = 100`.  No uniform proof yet.
        - FM-SEC54 (luna_max_jupiter): the same drop lemma and an exact
          angular form, with a proof on the subregion "every step at most
          `pi/2` and `mu_j <= 4 mu_i`".  On the centre line, `(HT*)` is
          `V(D_j - D_i)^2 >= D_j n_i^2`.  It also found the margin error
          corrected above.  FM-SEC59 continues on `(I)` with the
          local-metric chain.
        - FM-SEC59 (luna_max_jupiter): each local step of `(b)` is exactly
          `s_k = (D_k - D_(k+1)) / (2 sqrt(D_k (D_(k+1) +
          (X_k+1) A_(k+1)^2 / X_k^2)))`.  `(b)` fails on inner `q = 2`
          pairs (`(6,3,5,8)`, `(17,4,11,14)`, `(19,4,12,15)`,
          `(27,6,17,20)`), while `(b')` passes there.  At `(6,3,5,8)` the
          first local step turns by more than `pi/2` (det 1120,
          dot `-2016`).  `(I)` stays open at the shortest half-turns.
      - *`q = 2` plus sign* (`fm39/e_q2plus_split.py`,
        `fm39/e_q2plus_ratiofree.py`).  For all `a, e <= 60` (106,982
        pairs), every pair has `W >= 0` (plus sign trivial) or a positive
        definite plus-sign form (ratio-free); none has neither.  The
        condition `W >= 0` is not ratio-free: `W`'s own form is never PSD
        on the non-definite region.  So the actual ratio must be used
        there.
      - FM-SEC38 (luna_max_venus), Theorem OL's method on the `q = 2`
        plus sign (`fm39/e_q2plus_ol_repro.py`).
        - Proved: the terminal boundary and the definite region.
        - On the non-definite outer region, a continued-fraction ratio
          interval `n/d <= r <= R` plus three coefficient inequalities
          would finish; these are screened on the folded grid (904 outer
          cases pass), not proved.
        - The first non-definite case is `(a,e,j) = (19,3,20)`, where the
          form is indefinite but the actual value is 5614.  Open: all `a, e`
        (uniformity of the definite region in `(d, N, j)`, or ratio
        bounds as in Theorem OL's steps 2-5 for the non-definite pairs).
      - FM-SEC44 (luna_max_venus): outer-interval lemma.  With
        `n = N-j`, `v = n-1`, `h = 2j-N+2`, `U = j+2`, the gap-3 `W` is
        a binary quadratic `A + B r + C r^2` in `r = c_(j+1)/c_j`
        (explicit coefficients).  If `d^2 >= 4Uv` and
        `T = (h+1) d^2 - L - 3n(h-2)(h+v+3) >= 0`, then `W > 0` on
        Theorem OL's outer ratio interval `n/d <= r <= R`, so the plus
        sign is trivial there.  Proved.  The remaining step is the
        polynomial implication "non-definite plus form => `d^2 >= 4Uv`
        and `T >= 0`", screened on all 2,208 open non-definite cases with
        `a <= 60` (least `T = 1348` at `(25,3,24)`), not proved.
        Reproducer reruns exactly.
      - FM-SEC42 (luna_max_saturn): exact `t`-optimization of the single
        metric at `q = 2`.  `R(t) = F_j F_i - 4t(4V-t) Delta^2 =
        alpha t^2 + beta t + gamma`; if `alpha > 0` the optimum is
        `t* = -beta/(2 alpha)`.  Every `q = 2` pair with `W < 0` is closed
        by `t*` or an endpoint: 49,657 in the 600-row sample, and all of
        them in the grid `a, e <= 150` (1,681,577 `q = 2` pairs).
        Interior `t*` is needed, e.g. at `(a,e,j) = (5,181,112)`.  The
        70 earlier residuals all have `W > 0`.
      - FM-SEC50 (luna_max_venus), `(HT)` at `q = 2`: conditional.  The
        exact screen `3 <= e <= 150`, `e+2 <= a <= 150` gives 521,548
        pairs with `W < 0`, all closed by the interior `t*`.  All of them
        lie in the central regime `d^2 < 4(j+2)(N-j-1)`, where the
        outer-interval lemma does not apply.  So a proof must control the
        actual ratio there, plus a separate chart at `c_j = 0` (47 cases).
    - **Residual of (E) now** (`fm39/e_residual2.py`): 22,046 of 438,221
      pairs (5.03%), all with `q >= 2`: `q = 2` 8,364, `q = 3` 6,908,
      `q = 4` 3,948, `q = 5` 1,988, `q >= 6` 838.
    - **FM-SEC20 (luna_max_venus) on (ii):** proved at `e = 1` for every
      window, and for the reflection-plateau subcase at every odd `e`.
      - At `e = 1` `D` is log-concave, from the explicit form
        `D_k = C(a,k)^2 (a+1) Q_k/((k+1)(a-k+1)^2(a-k+2))` with
        `Q_k = (a-2k+1)^2 + a+1`.
      - The plateau subcase covers windows crossing the centre with
        `(N-2x+1) u >= 2x+C-N`.
      - The local log-concavity route fails, e.g. at `e = 3`, `a = 11`,
        `x = 8`, `C = 3`, where (ii) still holds.
      - Reproducer: `fm39/ii_e1_plateau_repro.py`.
    - **FM-SEC16 (luna_max_jupiter) on (ii):**
      - Proved: `C <= 1` (AM--GM), and centered windows `2x + C = N`.
        Centered windows follow from `D_(N-k) = D_k`, Theorem OL, and the
        plateau `D_(m-1) = D_m = D_(m+1)` at even `N`.
      - Exact screen: the 122,920 base windows reproduce; equality occurs
        only on centered `C = 2` plateaus with `a` odd.
      - Candidate Pascal step `a -> a+2`: with
        `D^(a+2)_(k+1) = D^(a)_k + L_k`, (ii) propagates if
        `sum_(k=x..x+C) L_k >= (C+1)(sqrt(D^(a+2)_(x+1) D^(a+2)_(x+C+1))
        - sqrt(D^(a)_x D^(a)_(x+C)))`.  Screened on 2,054,943 windows, no
        failure, unproved.
      - Killed (auxiliary criteria, no consumer use): log-concavity of the
        increments `L`, and of the ratios `D^(a+2)_(k+1)/D^(a)_k`.
      - *Consumer match.*  W at `C <= 2` is already proved by the
        certificates.  So (ii) is needed only for `C >= 3`, and the
        agent's first unresolved case `C = 2` is not consumed.
    - **FM-SEC13 (luna_max_venus) on (E):**
      - Proved: the `q = 0` boundary (`i = j+1`), where `W = D_j - D_(j+1)`.
        This is Theorem OL, plus a centre computation for `N = 2m` with
        `e` even:
        `D_m - D_(m+1) = t^2 (4(m+1)^2 - delta^2)/(2(m+1)^2(m+2))`.
      - Screen: (E) holds at all 438,221 pairs with `a, e <= 40`.  The
        least slack in the open region is 4, at `(a,e,j,i) = (3,5,7,9)`.
      - Kill (knob: using only positive semidefiniteness of `Gamma`, which
        the consumer does not restrict to): PSD alone does not give the
        Gram-path inequality, e.g. vectors `(1),(-1),(1)`.
      - No movement on `q >= 1`.
- **Four factors: split identity, closed form and Theorem LL4 (FM-SEC6,
  luna_max_venus; verified by the main agent, `mech/four_repro.py`,
  `mech/ll4_check.py`, `mech/four_explore.py`).**
  - *Split identity.*  `e = 2r-4` is even, so `W` is symmetric.  With
    shifted labels `L_i = kappa_i + 1`:
    `phi_r = tau(U_L1 U_L2 U_L3 U_L4) - sum_i W(prod_(j != i) U_Lj, U_Li)
    + sum_(3 pairings) W(U_Li U_Lj, U_Lk U_Ll)`.
    - The 2|2 pairing terms enter with a plus sign.
    - The main term is `>= 0`: `tau(U_l) = phi_(r-2)(hat S_l h_1^a)`, which
      is `>= 0` by Theorem OL for `r >= 3`.  For `r = 2` it is direct,
      from Catalan moments.
  - *Support and boundary.*  `W(U_p, U_q) = 0` for `p + q > N`, and
    `W(U_p, U_q) = c_q` for `p + q = N`.
  - **Theorem LL4.**  Sort `L_1 >= L_2 >= L_3 >= L_4` and put
    `mu_1 = L_4 + max(0, L_1 - L_2 - L_3)` and
    `mu_2 = L_1 - L_2 + L_3 - L_4`.  If `mu_1 > N` and `mu_2 > N`
    (`N = a + 2r - 4`), every cross term vanishes, and
    `phi_r(h_u h_v h_w h_z h_1^a) >= 0` for every `r >= 2`, `a >= 0`.
    It needs spread labels: four equal labels have `mu_2 = 0`.
  - *Closed form.*  `phi_r = M_4 - R_31 + R_22`, where
    `M_4 = sum_l (m_l - m_(l-2)) D_((N+l)/2)` and `R_31`, `R_22` are
    telescoped endpoint products `F_p(q) - F_p(q+2)`.
  - *Checks.*
    - The agent verifier reruns: split and closed form against direct
      Catalan values in 7,500 cases each, plus 1,792 kernel checks.
    - The main agent checked LL4 in 250 random cases with its own kernel
      evaluator, and the split identity in 141 large-label cases.
    - In those 141 cases the 2|2 pairing sum was negative in 47, but its
      magnitude was at most about `0.105 tau`.
  - Open: `M_4 - R_31 + R_22 >= 0` outside LL4, including equal or
    clustered labels.
- **Theorem LLm (spread labels, any number of factors; FM-SEC7,
  luna_max_venus; verified by the main agent, `mech/llm_check.py`).**
  - For `m <= 2r` general parts with shifted labels `L_i`, suffix `a`,
    `e = 2r - m` and `N = a + e`:
    `phi_r = sum over unordered splits {S, S^c} of (-1)^|S^c| W(U_S, U_(S^c))`,
    where `W(Q,P) = (-1)^e W(P,Q)`.
  - The main term `tau(prod U_Li) >= 0`, since
    - `tau(U_l) = phi_((e+1)/2)(h_(l-1) h_1^a)` (`e` odd), or
    - `tau(U_l) = phi_(e/2)(hat S_l h_1^a)` (`e >= 2` even), both `>= 0` by
      Theorem OL;
    - for `e = 0` it follows directly from Catalan moments.
  - Let `lambda(I)` be the least label `>= max(0, 2 max_I L - sum_I L)` of
    the parity of `sum_I L`; this is the smallest constituent of
    `prod_(i in I) U_(L_i)`.
  - **Theorem.**  If `lambda(S) + lambda(S^c) > N` for every proper split,
    every cross term vanishes by support, and `phi_r >= 0`.
  - It reduces to Theorem LL at `m = 3` (the support branch; the boundary
    branch adds the helpful `W(U_0, U_N) = -1`) and to Theorem LL4 at
    `m = 4`.
  - *Checks.*  The main agent checked it at `m = 5, 6` in 120 random
    cases with its own kernel evaluator, all nonnegative.  The agent's
    coverage count: 7.2% of the even-parity words on a small grid
    (`m <= 6`, parts `<= 4`, `r <= 4`, `a <= 4`).
  - Open: words with a split having `lambda(S) + lambda(S^c) <= N`, i.e.
    clustered labels.
- **Theorem LLm-S (spread labels, the whole consumer cone; FM-SEC8,
  luna_max_venus; verified by the main agent, `mech/llms_check.py`).**
  - Consumer words `w = prod_(i<=m) h_(kappa_i) prod_j hat S_(p_j) h_1^a` with
    `m <= 2r` and any number of `hat S`.
  - Absorb `(x-y)` into each `h`-factor, keep `hat S_p = U_p(x) + U_p(y)`,
    and expand over assignments of the factors to `x` or `y`:
    `phi_r = (1/2) sum_S (-1)^(|S^c cap H|) W(U_S, U_(S^c))`.
  - The one-sided terms give the main term `tau(prod U) >= 0`, by Theorem OL.
  - If `lambda(S) + lambda(S^c) > N` for every nontrivial assignment, every
    two-sided term vanishes by support, so `phi_r >= 0`.
  - It contains Theorem LLm.  At `r = 1`, FM53 and Lemma BP add cases it
    misses.
  - *Checks.*  The agent checked 5,580 words with one or two `hat S` against
    direct Catalan values.  The main agent checked 150 random words with
    1--4 `h` and 1--2 `hat S` factors satisfying the condition, with its own
    kernel evaluator; all were nonnegative.
- **Equal labels (FM-SEC9, luna_max_venus; checked by the main agent).**
  - *Squares.*  If `w = f^2` for a real word `f` and `a` is even, then
    `phi_r(w h_1^a) = (1/2) E[f^2 (x-y)^(2r) (x+y)^a] >= 0` pointwise.  If `a`
    is odd and `deg w` even, the value is 0 by sign reversal.  In
    particular `phi_r(h_k^m h_1^a) >= 0` for every even `m`.
  - *`k = 2`, every `m`.*  `h_2 = H_3 = x^2 + xy + y^2 - 2`, so the family
    is covered by the recorded FM18, FM20, FM22 and FM47.
  - `k = 1` is item (26) (all-ones), and `m = 1` is Theorem OL.
  - Open: odd `m >= 3` with `k >= 3`; the first case is
    `phi_2(h_3^3 h_1) = 3`.
- Direction-by-direction positivity fails.  For `(r,a,u,v,c) =
  (2,0,1,1,8)`, the `q = 0` radial integral is `-5/2431`, while the word
  is covered by Theorem THREE.  So only the full angular average can
  work.
- **FM-MECH2 open lead.**  The level step gives
  `phi_2(h_s h_t h_2 h_1^a) = [3 L_2(s,t,a+2) + L_3(s,t,a) - 8 L_2(s,t,a)]/4`,
  with `L_r(s,t,a) = phi_r(h_s h_t h_1^a)`.  With R3-TWO, the remaining
  sufficient inequality is `3 L_2(s,t,a+2) >= 8 L_2(s,t,a)`.  Not
  checked by the main agent.

**What this gives the full cone.**
- A precise, checked statement (qFM3_2) whose proof would close the whole
  `r = 2` H-only sector for every `kappa` at once.  Its remaining proof
  routes:
  - Kostka--Foulkes combinatorics (charge, cyclage);
  - a sign balance between the bigraded-pure blocks of the fusion Koszul
    complex (item (16)).
  - A vanishing theorem is ruled out: the fusion complex is not pure
    (item (10)), and slot complexes are not even (item (17)).
- It does not yet reach `r >= 3` or the `hat S` sectors.  The mechanism
  that makes integer levels positive at `r >= 3` is still unidentified.
- **Update (items (34)--(36)).**  Conjecture LP was the live target, and
  it is now KILLED (hooks `(14,1^14)`).  What survives, at integer `r`:
  - the closed form of FM3 (item (35));
  - the forced zeros;
  - stability, which reduces words with `d >= 1` to reweighted boundary
    words (item (36));
  - the local transport, which is feasible in every tested H-only case.
    An explicit transport rule is the current H-only target.
- **Superseded LP summary (kept for the record).**  The live full-cone
  target was Conjecture LP.
  - For every FM3 word, `phi_r(w)/phi_r(empty)` is an explicit rational
    function of `r`.
  - Its numerator is a product of forced integer zeros `1..d` (double,
    except `d`) and a factor `G_w`.
  - LP asks that `G_w > 0` on `[1, oo)`.  That one statement per word gives
    every level and every `hat S` sector at once.
  - Checked for 250+ words to total degree 10, and H-only to size 14.
  - It is tight on hooks `(k,1^k)`: the stray root tends to 1 from below.
