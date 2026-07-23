# Central-character Q3 search

## Candidate enlargement

For SU(2), let `chi_n` be the irreducible character of highest weight `n`.
All of these characters are real.  For SO(3), use the even-indexed characters
`chi_2n`.  The candidate cone is the uniform closure of the nonnegative span
of all these irreducible characters.  It strictly contains the cone generated
by the single adjoint character `chi_2`.

For a multiset `M` of minus-labelled characters and a multiset `P` of
plus-labelled characters, the tested quantity is

```text
J(M,P) = integral integral
         product_{i in M} (chi_i(g) - chi_i(h))
         product_{j in P} (chi_j(g) + chi_j(h)) dg dh.
```

Only even cardinalities of `M` need to be tested: odd cardinalities vanish
after interchanging `g` and `h`.

## Exact coefficient identity

Let `L` be the labelled list `M` followed by `P`.  For a subset `S` of its
positions, put

```text
N(S) = multiplicity of the trivial SU(2) module in tensor_{i in S} V_{L_i}.
```

Expanding the signed product and applying character orthogonality gives the
exact integer identity

```text
J(M,P) = sum_{S subset L} (-1)^{|S intersect M|} N(S) N(S^c).
```

The program `character_ring_iter/search_su2_central_q3.cpp` evaluates the
same integer a second way: it extracts the coefficient of `1 tensor 1` from

```text
product_{i in M} (chi_i tensor 1 - 1 tensor chi_i)
product_{j in P} (chi_j tensor 1 + 1 tensor chi_j)
```

in the representation ring of SU(2) x SU(2).  Both computations use the
Clebsch--Gordan rule

```text
chi_a chi_b = chi_|a-b| + chi_|a-b|+2 + ... + chi_a+b.
```

The executable checks the dynamic coefficient calculation against the direct
subset identity on 256 deterministic tuples before searching.

## A universal four-factor proposition

**Proposition 1.** Let `G` be a compact group and let every tested function be
a nontrivial real irreducible character of `G`.  Every Q3 signed-product
integral with at most four factors is nonnegative.

**Proof.** All-plus products are nonnegative after expansion because every
Haar integral of a product of characters is a trivial-representation
multiplicity.  An odd number of minus signs gives zero after interchanging
the two Haar variables.  It remains to consider two or four minus signs.

For two factors, Schur orthogonality gives `2 delta_ab`.  For three factors,
two minus signs and one plus sign give `2 N(a,b,c)`, where `N` is the
trivial multiplicity.  With four minus signs, direct expansion gives

```text
2 [N(a,b,c,d) + delta_ab delta_cd
                  + delta_ac delta_bd + delta_ad delta_bc],
```

which is nonnegative.  Finally, assign minus signs to `a,b` and plus signs
to `c,d`.  Expansion gives

```text
J/2 = N(a,b,c,d) + delta_ab delta_cd
                     - delta_ac delta_bd - delta_ad delta_bc.
```

If neither negative Kronecker product is one, the claim is immediate.  If
exactly one is one, the fourfold tensor product contains a trivial summand by
pairing the two equal `a/c` modules and the two equal `b/d` modules, so
`N(a,b,c,d) >= 1`.  If both are one, all four modules are equal; then the
positive Kronecker product is also one and again `N(a,b,c,d) >= 1`.  This
proves every remaining sign pattern.  QED.

Thus any counterexample made from real irreducible characters must use at
least five factors.  This proposition does not by itself cover a reducible
real generator formed as the sum of a complex character and its conjugate;
those generators are included in the SU(3) computation below.

In fact the support reduction below pushes this universal range one degree
farther.

**Proposition 1B (universal five-factor proposition).** Let `G` be a compact
group and let every tested generator be a nontrivial real irreducible
character.  Every Q3 signed-product integral with at most five factors is
nonnegative.

**Proof.** Proposition 1 handles at most four factors.  For five factors,
odd minus parity again gives zero.  If the same irreducible label occurs with
both signs, apply Proposition 5 to write the integral as a nonnegative linear
combination of signed integrals with one fewer factor; these are nonnegative
by Proposition 1.  It remains to consider disjoint plus and minus label
supports.

Let the five labelled positions be `L`.  Pair each subset `S` in the exact
coefficient identity with its complement.  Since the total minus parity is
even, the two terms have the same sign.  Choose from each complementary pair
the member of size at most two.  The empty-set contribution is `N(L)`.  Every
singleton contribution vanishes because its label is a nontrivial
irreducible.  A two-element contribution is

```text
(-1)^(|S intersect M|) delta_(L_i,L_j) N(L minus S).
```

It is nonnegative when the two positions have the same sign.  When their
signs differ, disjoint support makes the Kronecker delta zero.  Hence the
whole integral, which is twice the sum of these representatives, is
nonnegative.  QED.

Thus a counterexample in this setting would require at least six factors.
At six factors the complementary middle layer consists of `3|3` splits;
their invariant multiplicities can be nonzero without any repeated label,
so the preceding universal argument genuinely stops there.  For `SU(2)`,
however, the low-spin channels dominate this middle layer as follows.

**Proposition 1C (the complete six-factor SU(2) sector).** Every Q3 signed
product of at most six nontrivial irreducible `SU(2)` characters is
nonnegative.

**Proof.** Only six factors require consideration.  Apply Proposition 5 and
induction on the number of factors until the plus and minus label supports
are disjoint.  Pair complementary subsets in the exact coefficient identity.
The empty-set term contributes `N(L)`, singleton terms vanish, and every
two-element term is nonnegative: opposite-sign equalities are excluded by
disjoint support, while same-sign equalities have positive sign.  The only
remaining terms are the ten unordered `3|3` splits.  A threefold `SU(2)`
invariant multiplicity is either zero or one.  Consequently it is enough to
prove

```text
N(p_1,...,p_6) >= T(p_1,...,p_6),                (P1C.1)
```

where `T` is the number of unordered splits into two invariant-admissible
triples.

For a triple `(a,b,c)`, let `m_r(a,b,c)` be the multiplicity of `V_r` in
`V_a tensor V_b tensor V_c`.  If the triple is invariant-admissible, order
it so that `c` is its largest label.  With

```text
I(x,y)={|x-y|,|x-y|+2,...,x+y},
```

the Clebsch--Gordan rule gives

```text
m_r(a,b,c)=|I(a,b) intersect I(c,r)|.
```

Admissibility says `c in I(a,b)`.  Since all three labels are positive,
`I(a,b)` has at least two entries.  The neighboring-entry description of the
two intervals therefore gives

```text
m_0 >= 1,       m_2 >= 2,       m_4 >= 1.        (P1C.2)
```

If all three labels are even, they are at least two and `I(a,b)` has at
least three entries.  The same endpoint check, including the sole smallest
case `(2,2,2)`, gives the stronger bounds

```text
m_0 >= 1,       m_2 >= 2,       m_4 >= 2,
m_6 >= 1.                                      (P1C.3)
```

Choose any admissible split `A|B` and associate the sixfold product according
to that split.  Orthogonality gives

```text
N(p_1,...,p_6)=sum_(r>=0) m_r(A)m_r(B).
```

Thus `(P1C.2)` implies `N>=1+4+1=6`.  If all six labels are even,
`(P1C.3)` implies `N>=1+4+4+1=10`.

It remains only to compare these bounds with `T`.  Let `o` be the number of
odd labels.  Each admissible triple contains an even number of odd labels.
If `o=2`, there are at most four parity-compatible `3|3` splits; if `o=4`,
there are at most six; and if `o=6`, there are none.  Hence whenever at least
one label is odd, `T<=6<=N`.  If `o=0`, all labels are even and there are only
ten unordered `3|3` splits in total, so `T<=10<=N`.  This proves `(P1C.1)`.

Returning to the signed subset expansion, discard the nonnegative pair and
positive-sign middle terms and bound the number of negative middle terms by
`T`.  The remaining half-integral is at least `N-T>=0`; doubling proves the
claim.  QED.

The exact C++ diagnostic
`character_ring_iter/search_su2_six_factor_middle.cpp` independently checks
the six-factor subset decomposition, `(P1C.1)`, and the stronger intermediate
counts.  With labels at most twelve it passed all 130,780 disjoint-support
signed cases; the all-minus portion covers all 12,376 sorted label sextuples.

Thus any `SU(2)` counterexample must use at least seven factors.

## The seven-factor straightening target

For seven factors, pair complementary subsets and retain one representative
of every `3|4` split.  After the disjoint-support reduction, every nonzero
pair term has equal labels of the same sign and is nonnegative.  If `N_7` is
the sevenfold invariant multiplicity, `P` is the sum of these equal-pair
terms, and `T_-` and `T_+` are the negative- and positive-sign `3|4`
contributions, respectively, the exact half-integral is

```text
N_7+P+T_+-T_-.
```

Exact enumeration through label ten, 138,390 disjoint signed cases, gives

```text
T_- <= N_7+P.                                  (7.1)
```

Raw domination by `N_7` is false: minus labels `[2,2,2,3]` and plus labels
`[1,1,1]` give `N_7=17` and `T_-=19`.  The pair term is `24`, so `(7.1)`
still has substantial slack.  A stronger bound by all middle terms also
fails, for example at six minus labels `2` and one plus label `4`.

There is now a precise standard-monomial explanation for the correction
`P`.  Realize `V_p=Sym^p(C^2)`.  The `SL_2` invariant ring is generated by
brackets `[ij]`, with quadratic Plucker relations, and its noncrossing
bracket monomials form a basis.  In multidegree `(p_1,...,p_7)` this basis
has cardinality `N_7`.  An invariant-admissible triple has a unique triangle
bracket monomial.  Multiplying it by each noncrossing invariant monomial on
the complementary four vertices gives exactly the objects counted by its
`3|4` contribution.  Thus the `T_-` negative objects straighten into the
`N_7`-dimensional global invariant space.

Let `rho` be the rank of these straightened negative products and put

```text
kappa=T_--rho.
```

The seven-factor inequality follows from the single syzygy bound

```text
kappa <= P,                                    (7.2)
```

because `rho<=N_7`.  This is sharper structurally than `(7.1)`: it says that
only relations among the factorized `3|4` products need boundary
compensation.  When all seven labels are distinct, `P=0`, so `(7.2)` asserts
linear independence of the negative products.  Repeated equal labels supply
the only proposed relation generators, through a full pair bracket times an
invariant on the remaining five vertices.

The strict C++ program
`character_ring_iter/search_su2_seven_factor_straightening.cpp` constructs
all degree-constrained noncrossing monomials, performs Plucker straightening,
and computes the source rank modulo `1,000,000,007`.  A modular rank is a
lower bound for the rational rank, so a verified modular kernel bound proves
the corresponding rational bound.  With labels at most seven it checked all
12,691 disjoint signed cases and found `kappa<=P` in every case.  It included
active distinct-label cases, where the kernel was zero.  In the first raw-
domination counterexample above, the 19 negative products have modular rank
16 and kernel dimension three, well below `P=24`.  This is an exact finite
certificate and identifies `(7.2)` as the next algebraic lemma; it is not yet
a proof for unbounded labels.

A simple strengthening via one trivalent-tree valuation is false.  Among all
compatible four-split tree valuations, the case with minus labels `[3,5]`
and plus labels `[1,2,4,6,7]` has at least two colliding source valuations.
Here `P=0`, but the exact modular rank is full (`22` sources, kernel zero), so
the failure concerns only the one-tree separation shortcut, not `(7.2)`.
Any valuation proof must combine charts or retain lower terms; alternatively
one can prove exactness of the proposed Plucker-syzygy complex directly.

Exact nullspace extraction sharpens this warning.  For minus labels `[1,1]`
and plus labels `[2,2,3,3,4]`, the 24 source products have rank 19.  A modular
nullspace basis has five relations, all with coefficients in `{0,1,-1}` and
all alternating under equal-label transpositions.  Relations are not confined
to repeated fundamentals: with six minus labels `2` and one plus label `4`,
the 40 sources have rank 16 and kernel dimension 24.  The pair term is 90.
Thus equal-label symmetry remains the correct locus, but a proof must handle
all repeated irreducibles, not only label `1`.

Two natural ways to turn that observation into `(7.2)` are false.  First,
projecting the half of each transposition pair to the singlet channel does not
separate the kernel: its rank on the preceding five-dimensional kernel is only
two (and is 19 on the 24-dimensional all-label-`2` kernel).  Second, changing
the circular order of the seven vertices does not restore a leading-monomial
proof.  In the distinct-label case with minus `[3,5]` and plus
`[1,2,4,6,7]`, all 720 circular orders fixing one vertex have at least two
leading collisions even though the exact source rank is full.  Any successful
standard-monomial argument therefore needs a genuine multi-chart or
lower-term determinant calculation.

There is also one complete infinite ray that follows directly from parity.

**Proposition 1A (every odd-character ray).** Fix an odd label `p`.  Every
signed-product integral made only from the `SU(2)` character `chi_p` is
nonnegative, with no bound on the number of factors.

**Proof.** Write `x=chi_p(g)` and `y=chi_p(h)`.  With `2m` minus factors and
`q` plus factors the integrand is

```text
(x-y)^(2m) (x+y)^q.
```

It is pointwise nonnegative when `q` is even.  Since `p` is odd, simultaneous
central translation by `-I` sends `(x,y)` to `(-x,-y)`.  It preserves Haar
measure and negates the integrand when `q` is odd, so that integral is zero.
QED.

This proves full `GKS2*` on every one-generator odd-character ray, including
the fundamental ray.  It does not extend by the same pointwise argument to
even-indexed characters, whose Haar distributions are not symmetric about
zero, or to mixtures of character labels.

## Executed domains

The following exact runs were completed on 2026-07-21.

1. The persistent C++ program exhaustively tested every separately sorted
   minus and plus multiset with at most eight total factors and labels
   `1,...,6` for SU(2), or `2,4,...,12` for SO(3).  This is 60,087 cases per
   group.  No value was negative.
2. A direct Python implementation of the subset identity exhaustively tested
   1,510,795 SU(2) cases with labels `1,...,10` and at most eight factors.
   It separately tested the same 60,087-case SO(3) domain.  No value was
   negative.
3. The C++ program then tested 50,000 deterministic invariant-admissible
   tuples for each group, with 4--18 factors.  SU(2) weights were at most 20
   and SO(3) weights at most 40.  The minimum in each 50,000-case sample was
   exactly 2.
4. An independent SU(3) Weyl-constant-term search tested 5,975 signed
   products from the four self-dual irreducibles `chi_(p,p)`, `1 <= p <= 4`,
   and 20,655 products from five real generators consisting of the adjoint
   character and four conjugate-pair sums.  Every tuple had at most eight
   factors.  No value was negative.  The implementation is
   `character_ring_iter/search_su3_real_central_q3.py`; before searching, it
   reproduces the ten-term SU(3) adjoint moment prefix used by the main paper.

The admissibility filter retains only tuples whose total tensor product can
contain the trivial SU(2) module: the total highest weight is even and the
largest highest weight is at most the sum of the others.  A tuple failing
this condition has identically zero `J`, so excluding it cannot hide a
negative value.

To rebuild and run the persistent search:

Set `THREADS` to the largest value compatible with the available RAM and the
current machine load.  On a dedicated machine with adequate memory,
`THREADS=$(nproc)` uses all available hardware threads; reduce it when other
jobs need capacity or per-thread memory becomes limiting.

```sh
g++ -O3 -DNDEBUG -std=c++20 -fopenmp \
  -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Werror \
  character_ring_iter/search_su2_central_q3.cpp \
  -o /tmp/search_su2_central_q3
THREADS="${THREADS:-$(nproc)}"
OMP_NUM_THREADS="$THREADS" /tmp/search_su2_central_q3 6 8 20 18 50000
python3 character_ring_iter/search_su3_real_central_q3.py
```

## What these calculations do and do not establish

They rule out a large class of small and medium central-character
counterexamples.  They do not prove Q3 for the infinite central character
cone.

There is relevant but insufficient adjacent theory.  Gasper proves
nonnegative linearization coefficients for the applicable Jacobi and
ultraspherical polynomials in *Linearization of the Product of Jacobi
Polynomials II*, Theorem 1
(<https://doi.org/10.4153/CJM-1970-065-4>).  Bakry and Huet explain that
character bases have both the GKS and hypergroup properties and formulate a
GKS2 conjecture from those properties in Section 2.6 of *The Hypergroup
Property and Representation of Markov Kernels*
(<https://arxiv.org/abs/math/0601605>).  Nonnegative product linearization,
the hypergroup property, and GKS2 are not the full signed-product condition
Q3, so neither source supplies the missing theorem.

## Current proof target

For every finite list of SU(2) highest weights, define the symmetric set
function

```text
w(S) = N(S) N(S^c).
```

The rank-one central-character theorem is equivalent to proving that every
even Walsh coefficient

```text
sum_S (-1)^{|S intersect M|} w(S)
```

is nonnegative.  The next proof step is therefore a nonnegative path or
Temperley--Lieb expansion of these Walsh coefficients.  Without such an
identity, the infinite-cone statement remains a conjecture.

## An exact rank-one reduction

Write `x = chi_1` and use the convention `chi_{-1} = 0`.  The SU(2)
recurrence

```text
x chi_j(x) = chi_{j+1}(x) + chi_{j-1}(x)
```

gives the telescoping Christoffel--Darboux identity

```text
chi_n(x) - chi_n(y)
  = (x-y) D_n(x,y),

D_n(x,y)
  = sum_{j=0}^{n-1} chi_j(x) chi_{n-1-j}(y).
```

In particular, `D_n` has nonnegative coefficients in the product-character
basis.  If there are `2m` minus-labelled factors, the full integrand is

```text
(x-y)^(2m)
  product_{i in M} D_i(x,y)
  product_{j in P} (chi_j(x) + chi_j(y)).
```

Everything after `(x-y)^(2m)` has nonnegative product-character
coefficients.  Thus the only remaining obstruction is the common Weyl
difference factor.

That obstruction has an exact parity description.  If `mu` is the SU(2)
class measure, then

```text
K_m(a,b)
  = integral integral (x-y)^(2m) chi_a(x) chi_b(y) dmu(x)dmu(y)
  = (-1)^b L_m(a,b),

L_m(a,b)
  = integral integral (x+y)^(2m) chi_a(x) chi_b(y) dmu(x)dmu(y) >= 0.
```

The equality follows by `y -> -y` and
`chi_b(-y) = (-1)^b chi_b(y)`.  The last inequality follows by expanding
`(x+y)^(2m)`: every one-variable integral
`integral x^r chi_a(x) dmu(x)` is the multiplicity of the trivial module in
`V_1^(tensor r) tensor V_a`.  Consequently, if the nonnegative remainder is

```text
sum_{a,b} c_{a,b} chi_a(x) chi_b(y),
```

then Q3 is exactly the assertion that its even-`b` contribution dominates
its odd-`b` contribution after weighting by `L_m(a,b)`.  Positivity of the
Christoffel--Darboux factors alone does not prove this last dominance.

## Two coefficientwise positivity routes that fail

The all-minus remainder has the symmetric-function form

```text
product_i D_(n_i)(x,y) = product_i h_(n_i-1)(z,z^-1,w,w^-1).
```

Termwise positivity after expanding this product in Schur functions would be
sufficient, but the required functional is not Schur-positive.  The exact C++
calculation in `character_ring_iter/search_su2_schur_kernel.cpp` gives

```text
integral integral (x-y)^8 s_(1,1)(z,z^-1,w,w^-1) dmu(x)dmu(y) = -168.
```

Here `s_(1,1)=2+xy`, and direct expansion using the semicircle Catalan moments
independently gives the same integer.

A second possible proof would split an all-minus product into two odd-sized
halves and expand each half in the antisymmetric orthogonal basis

```text
chi_p(x)chi_q(y) - chi_q(x)chi_p(y),  p > q.
```

There is no label-independent orthant containing all such coefficient
vectors.  For the `(p,q)=(2,1)` coefficient, the exact values are `-3` for
the label triple `[1,1,1]` and `+1` for `[1,2,2]`.  The OpenMP C++ search
`character_ring_iter/search_su2_half_product_cone.cpp` found 1,195 failures
among the 1,540 sorted triples with labels at most 20.  Thus any Gram-matrix
argument must group several recoupling channels rather than prove positivity
coefficient by coefficient.

## Why the first diagram argument fails

For fundamental factors, `N(S)` counts noncrossing matchings on `S`.  After
overlaying a matching on `S` with one on `S^c`, forgetting the two colors
gives a chord matching whose crossing graph must be bipartite.  It is
tempting to sum over the two colorings of each crossing component and claim
that every such orbit is nonnegative.  This is false in the first possible
case.

Take four fundamental factors, put minus signs at vertices `1,2`, and use
the crossing matching `(1,3),(2,4)`.  Its two valid edge colorings each put
exactly one minus vertex in the first color, so each has sign `-1`; the orbit
contributes `-2`.  The complete signed integral is nevertheless `+2`.
Therefore any diagram proof must include an uncrossing/skein cancellation
between different uncolored matchings; color switching within one matching
cannot prove the theorem.

There is a related operator reformulation, but its strongest tempting form is
also false.  Put

```text
H = Inv(V_(L_1) tensor ... tensor V_(L_n)).
```

For each subset `S`, the subspace

```text
Inv(tensor_(i in S) V_(L_i))
  tensor Inv(tensor_(i notin S) V_(L_i))
```

embeds canonically in `H`; let `P_S` be its orthogonal projector.  Then

```text
rank(P_S) = N(S)N(S^c),
J(M,P) = trace_H A_M,
A_M = sum_S (-1)^(|S intersect M|) P_S.           (OP)
```

Thus positivity of `A_M` would prove the full theorem.  It does not hold.
For four fundamental factors and minus positions `{1,2}`, let `Q_12|34`,
`Q_13|24`, and `Q_14|23` be the rank-one projectors onto the three normalized
pair-singlet vectors.  On the two-dimensional global singlet space,

```text
A_M = 2[I + Q_12|34 - Q_13|24 - Q_14|23].
```

The three pair-singlet lines have mutual overlap of absolute value `1/2`, so
their projectors sum to `(3/2)I`.  Consequently

```text
A_M = 4 Q_12|34 - I,
spectrum(A_M) = {3,-1},
trace(A_M) = 2.
```

The required trace is positive even though the operator is indefinite.  The
C++ program `character_ring_iter/test_su2_projection_psd.cpp` independently
constructs the singlet projectors in the spin basis and numerically reproduces
the eigenvalues.  Hence a proof may use the exact trace identity `(OP)`, but
not the stronger positive-semidefinite claim.

## Relation to the classical GKS2* problem

Bakry and Echerbault call precisely the condition that every signed monomial
in the basis sums and differences have nonnegative doubled integral the
`GKS2*` property.  Thus the proposed all-character SU(2) enlargement is an
instance of their full `GKS2*` problem, not a consequence of nonnegative
linearization.  They prove `GKS2*` for abelian character systems and several
small finite-group systems, but state that they have no general mechanism:
*Sur les inegalites GKS*, Section 3
(<https://www.numdam.org/item/SPS_1996__30__178_0/>).  Bakry and Huet later
state the weaker implication `GKS + HGP => GKS2` as a conjecture:
*The Hypergroup Property and Representation of Markov Kernels*, Section 2.6
(<https://arxiv.org/abs/math/0601605>).  This explains why the final parity
dominance above is the substantive step.

As a finite approximation check, the exact `SU(2)_k` fusion-rule calculation
found no negative value for levels `k=2,...,10`, all nontrivial labels, and at
most eight factors.  The respective exhaustive case counts were `210`,
`1,354`, `5,975`, `20,655`, `60,087`, `153,615`, `355,113`, `757,185`, and
`1,510,795`.  These finite systems support, but do not imply, the ordinary
SU(2) statement.  The reproducible exact C++ implementation is
`character_ring_iter/search_su2_level_central_q3.cpp`; it independently
checks 128 cases against direct subset expansion before searching.  A strict
build and the extended run are:

```sh
g++ -O3 -DNDEBUG -std=c++20 -fopenmp \
  -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Werror \
  character_ring_iter/search_su2_level_central_q3.cpp \
  -o /tmp/search_su2_level_central_q3
THREADS="${THREADS:-$(nproc)}"
OMP_NUM_THREADS="$THREADS" /tmp/search_su2_level_central_q3 10 8
```

## Exact stabilization of the fusion approximation

**Proposition 2 (pathwise stabilization).** Let `a_1,...,a_r` be ordinary
SU(2) highest weights.  Suppose that every `a_i <= k` and

```text
a_1 + ... + a_r <= 2k.
```

Then the multiplicity of the tensor unit in the ordinary tensor product is
equal to its multiplicity in the `SU(2)_k` fusion product.

**Proof.** The ordinary multiplicity counts paths
`c_0=0,c_1,...,c_r=0` for which the triple
`(c_(i-1),a_i,c_i)` satisfies the triangle inequalities and has even sum.
The level-`k` fusion rule imposes the additional inequality

```text
c_(i-1) + a_i + c_i <= 2k.
```

Along any ordinary path ending at zero, the triangle inequalities applied
from the two ends give

```text
c_(i-1) <= a_1 + ... + a_(i-1),
c_i     <= a_(i+1) + ... + a_r.
```

Adding these bounds and `a_i` proves the additional fusion inequality.
Moreover, `c_i` is at most both the prefix sum through `i` and the remaining
suffix sum, so `c_i <= (a_1+...+a_r)/2 <= k`.  Thus every ordinary path is a
level-`k` path.  The reverse inclusion follows directly because the fusion
rule only truncates the ordinary Clebsch--Gordan outputs.  QED.

**Corollary 3 (signed-integral stabilization).** For a fixed signed list of
weights with total weight `T` and maximum weight `A`, its ordinary Q3 integer
is exactly its `SU(2)_k` fusion Q3 integer for every

```text
k >= max(A, ceil(T/2)).
```

**Proof.** In the subset formula, both `N(S)` and `N(S^c)` involve total
weight at most `T`, and all of their labels are at most `A`.  Proposition 2
therefore applies to every summand.  QED.

Consequently, a proof of `GKS2*` for every finite system `SU(2)_k`, uniform
only in the sense that it covers every level, would imply the ordinary
SU(2) signed-character statement without a limiting argument: each fixed
integer stabilizes at the displayed finite level.  Uniform closure of the
nonnegative character cone is then handled separately by continuity of each
fixed signed-product integral in the sup norm.

The stabilization claim was also checked independently by
`character_ring_iter/verify_su2_fusion_stabilization.cpp`.  A strict C++ build
tested all 60,087 separately sorted signed lists with labels at most six and
at most eight factors, both at the displayed threshold and one level above;
every fusion integer equalled the ordinary integer:

```sh
g++ -O3 -DNDEBUG -std=c++20 -fopenmp \
  -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Werror \
  character_ring_iter/verify_su2_fusion_stabilization.cpp \
  -o /tmp/verify_su2_fusion_stabilization
THREADS="${THREADS:-$(nproc)}"
OMP_NUM_THREADS="$THREADS" /tmp/verify_su2_fusion_stabilization 6 8
```

## Boundary domination for the two-minus sector

**Proposition 4 (exact boundary formula).** Fix a level `k`, with the ordinary
SU(2) ring interpreted as `k=infinity`.  Expand a product of plus factors as

```text
product_j (e_(p_j) tensor 1 + 1 tensor e_(p_j))
  = sum_(u,v) c_(u,v) e_u tensor e_v.
```

For two minus labels `a,b`, let `a star_k b` denote the (possibly truncated)
Clebsch--Gordan interval.  Then the signed Q3 coefficient is

```text
J_k(a,b;P)
  = sum_(r in a star_k b) (c_(r,0) + c_(0,r))
      - c_(a,b) - c_(b,a)
  = 2 (sum_(r in a star_k b) c_(r,0) - c_(a,b)).
```

**Proof.** Multiply the displayed plus product by

```text
(e_a tensor 1 - 1 tensor e_a)
(e_b tensor 1 - 1 tensor e_b).
```

The two same-side terms contribute the two boundary sums because the tensor
unit occurs in `e_r e_s` exactly when `r=s`.  The two opposite-side terms
contribute `-c_(a,b)-c_(b,a)`.  The coefficient array is symmetric because it
starts at `(0,0)` and every plus factor is invariant under exchanging the two
coordinates, giving the second equality.  QED.

Thus the complete two-minus sector at level `k` is equivalent to the named
boundary-domination estimate

```text
BD_k(a,b;P):  c_(a,b) <= sum_(r in a star_k b) c_(r,0).
```

The function `q3_two_minus_boundary` in
`character_ring_iter/search_su2_level_central_q3.cpp` computes the right-hand
side independently of the signed-product update.  The executable checks the
formula on 128 deterministic finite-level cases in addition to its 128 direct
subset cross-checks.

Additional exact random runs at level `50`, completed on 2026-07-21, tested
100,000 deterministic samples with at most 24 factors.  The final run's
minimum was zero, attained for minus labels `[1,1]` and plus labels
`[1,2,2,7,12,16,17,17,21,22,31,37,39,45,49,49]`; no tested value was
negative.  This is falsification evidence, not a proof of boundary domination
or of full `GKS2*`.

## Two exact reductions of the remaining target

**Proposition 5 (disjoint-support reduction).** In either the ordinary
`SU(2)` character ring or any finite fusion ring `SU(2)_k`, it is enough to
verify `GKS2*` for signed monomials in which no nontrivial character label
occurs among both the plus factors and the minus factors.

**Proof.** Write

```text
s_i(x,y) = chi_i(x) + chi_i(y),
d_i(x,y) = chi_i(x) - chi_i(y).
```

If a signed monomial contains both `s_i` and `d_i`, remove one copy of each
and use

```text
s_i d_i
  = chi_i(x)^2 - chi_i(y)^2
  = sum_r N_(i,i)^r d_r.
```

All ordinary tensor-product coefficients, and all level-`k` fusion
coefficients, in the last sum are nonnegative.  Every resulting monomial has
one fewer factor.  Repeating terminates and expresses the original integral
as a nonnegative linear combination of integrals with disjoint plus/minus
label supports.  The term `r=0` vanishes because `d_0=0`.  QED.

This is the support reduction in Bakry--Echerbault, Remark 3.7, but the
displayed fusion-ring calculation proves the precise version needed here
without importing a finite-space normalization convention.

There is a second exact reduction specific to the `SU(2)` parity grading.

**Proposition 5A (odd-label sign toggle).**  In the ordinary ring or any
finite fusion ring `SU(2)_k`, simultaneously interchange plus and minus signs
on every odd-labelled factor.  The signed contraction is unchanged.  Thus,
if

```text
M'=M_even multiset-union P_odd,
P'=P_even multiset-union M_odd,
```

then

```text
J(M,P)=J(M',P').                                      (P5A.1)
```

**Proof.**  The map

```text
J_2(e_a tensor e_b)=(-1)^b e_a tensor e_b
```

is induced by central translation by `-I` in the second `SU(2)` variable.
It is a ring automorphism and preserves the trivial-coefficient functional;
the same statements hold in `SU(2)_k` because fusion respects label parity.
For even `i`, both `s_i` and `d_i` are fixed.  For odd `i`, one has

```text
J_2(s_i)=d_i,                 J_2(d_i)=s_i.
```

Applying `J_2` to the entire signed product and then taking its trivial
coefficient proves `(P5A.1)`.  QED.

Consequently sign patterns need only be studied modulo complementation on
the odd-labelled positions.  In particular, one may choose the representative
with no more odd minus positions than odd plus positions.  If the total
number of odd-labelled factors is odd, simultaneous central translation in
both variables already makes the contraction zero.

For the second reduction, let `N_p` be the symmetric fusion matrix of label
`p`, let `e_p` be the corresponding coordinate vector, and define the fusion
Hankel extension

```text
H(v) = sum_r v_r N_r.
```

If `C_P=(c_(a,b))` is the coefficient matrix of a plus product, put

```text
d_P = C_P e_0,
E_P = H(d_P) - C_P.
```

By Proposition 4, simultaneous boundary domination for every pair `(a,b)`
is exactly entrywise nonnegativity of `E_P`.

**Proposition 6 (closed defect dynamics).** Starting with the empty plus
product, the boundary defect satisfies

```text
E_empty = I - e_0 e_0^T,
E_(pP) = L_p(E_P),
L_p(E) = N_p E + E N_p - H(E e_p).
```

Every `E_P` is symmetric and has zero row and column at label zero.

**Proof.** The empty coefficient matrix is `e_0 e_0^T`, whose boundary is
`e_0`; since `H(e_0)=N_0=I`, the initial formula follows.  Adding a plus
factor gives

```text
C_(pP) = N_p C_P + C_P N_p,
d_(pP) = N_p d_P + C_P e_p.
```

Fusion associativity gives `H(N_p v)=N_p H(v)`.  Also
`H(d_P)e_p=N_p d_P`.  Substituting `C_P=H(d_P)-E_P` into the last two
displays and cancelling the two copies of `N_p H(d_P)` gives exactly

```text
H(d_(pP)) - C_(pP)
  = N_p E_P + E_P N_p - H(E_P e_p).
```

Symmetry is preserved.  Multiplying the displayed recurrence by `e_0` and
using `N_p e_0=e_p` and `H(v)e_0=v` shows `E_(pP)e_0=0`; the initial defect
has the same property.  QED.

Equivalently, on symmetric matrices the map

```text
delta(C) = H(C e_0) - C
```

intertwines the plus update with `L_p`.  Its kernel is precisely the
`(k+1)`-dimensional fusion-Hankel subspace `H(v)`.  Hence the defect space has
dimension

```text
(k+1)(k+2)/2 - (k+1) = k(k+1)/2.
```

The remaining two-minus theorem is now the concrete orbit statement

```text
L_(p_n) ... L_(p_1)(I-e_0 e_0^T) >= 0 entrywise
```

for every level and every word of plus labels (and, by Proposition 5, it is
enough to pair the resulting entries with minus labels absent from that
word).  The subtraction in `L_p` shows why the full nonnegative matrix cone
is not invariant.  A proof therefore requires a smaller invariant cone or a
positive basis for this particular defect orbit.

## A two-layer switching reduction in the ordinary ring

The ordinary Clebsch--Gordan intervals give a sharper reduction that respects
highest-weight parity.  For `a,b >= 2`, define

```text
Q_P(a,b) = E_P(a,b) - E_P(a-2,b-2).
```

Since the interval for `(a-2) star (b-2)` is obtained from the interval for
`a star b` by deleting its two largest elements, one has the exact formula

```text
Q_P(a,b)
  = d_P(a+b) + d_P(a+b-2)
      - c_P(a,b) + c_P(a-2,b-2).                 (TLS)
```

At the parity-chain edge, set

```text
Q_P(a,1) = E_P(a,1)
         = d_P(a-1) + d_P(a+1) - c_P(a,1).
```

The proposed smaller inequality is the following.

```text
Two-layer switching (ordinary SU(2)):
If neither a nor b occurs among the plus labels P, then
Q_P(a,b) >= 0,
with the displayed edge convention when min(a,b)=1.
```

**Proposition 7 (two-layer switching implies the two-minus sector).** If the
two-layer switching inequality holds for every ordinary plus word, then every
ordinary `SU(2)` signed product with exactly two minus factors is
nonnegative.

**Proof.** Induct on the number of plus factors.  The empty plus word has
`E=I-e_0e_0^T`, so the assertion holds.  For a nonempty word `P`, first
suppose that a minus label, say `a`, occurs in `P`.  Pair one corresponding
plus factor with that minus factor and apply Proposition 5:

```text
s_a d_a = sum_r N_(a,a)^r d_r.
```

The desired two-minus integral becomes a nonnegative sum of two-minus
integrals with one fewer plus factor, and these are nonnegative by the outer
induction hypothesis.

It remains to take `a,b` absent from the support of `P`.  If one label is
zero the defect is zero, and if the smaller label is one the edge inequality
gives `E_P(a,b)>=0`.  Otherwise `(TLS)` gives

```text
E_P(a,b) >= E_P(a-2,b-2).
```

Repeat down the parity diagonal.  If the chain reaches label zero or one,
the preceding edge argument applies.  If it first reaches a pair containing
a label from the support of `P`, Proposition 5 and the outer induction
hypothesis show that this terminal defect is nonnegative.  Adding the
nonnegative two-layer increments proves `E_P(a,b)>=0`.  Proposition 4 then
gives the signed integral.  QED.

This reduction uses the representation-ring identity

```text
e_a e_b = e_(a-2)e_(b-2) + e_(a+b) + e_(a+b-2)
```

for `a,b>=2`: the last two summands are precisely the two outer
Clebsch--Gordan layers.  Thus the remaining two-minus problem is a switching
injection into the lower pair `(a-2,b-2)` and those two outer channels.

There is now a second, inductive reduction of the interior part of `(TLS)`.
For a symmetric defect matrix `E`, a prospective next plus label `p`, and
`a>=b>=2`, define the two one-sided increments

```text
A_p^E(a,b)=(N_p E)_(a,b)-(N_p E)_(a-2,b-2),
B_p^E(a,b)=(E N_p)_(a,b)-(E N_p)_(a-2,b-2).
```

For a nondecreasing prefix word `P`, write `E=E_P` and take
`p>=max(P)`.  The new exact target is

```text
A_p^E(a,b) >= E_(a+b,p),
B_p^E(a,b) >= E_(a+b-2,p),                      (OSD)
```

whenever `a,b` are absent from the support of `P` and differ from `p`.

**Proposition 7A (one-sided domination closes the interior induction).**
Assume `(OSD)` for every ordered prefix.  If the prefix defect is entrywise
nonnegative, then adjoining its next largest plus label preserves `(TLS)` at
every disjoint target `a>=b>=2`.

**Proof.** Proposition 6 gives

```text
E_(pP)=N_p E+E N_p-H(E e_p).
```

Take the difference at `(a,b)` and `(a-2,b-2)`.  The first two terms are
exactly `A_p^E(a,b)+B_p^E(a,b)`.  For every vector `v`, the difference

```text
H(v)_(a,b)-H(v)_(a-2,b-2)
```

is `v_(a+b)+v_(a+b-2)`, because the lower Clebsch--Gordan interval is
obtained by deleting the two outer channels.  With `v=E e_p`, this gives the
exact identity

```text
Q_(pP)(a,b)
 = A_p^E(a,b)+B_p^E(a,b)
     -E_(a+b,p)-E_(a+b-2,p).                   (7A.1)
```

The two inequalities `(OSD)` make the right side nonnegative.  QED.

This formulation automatically includes the reflected cases `a<p` or
`b<p`; no informal nesting of Clebsch--Gordan intervals is being assumed.
In the stable region the one-sided increments become sums of old `Q`
increments, while below `p` the exact reflected boundary entries provide the
additional terms.  The strict C++ mode
`character_ring_iter/search_su2_defect_cone.cpp --local` checked `(OSD)` for
all 1,716 nondecreasing prefix words with labels at most six and at most seven
factors, every next label `p` between the largest prefix label and six, and
every disjoint interior target.  All cases passed.  This is a sharper
inductive target than `(TLS)`, but it remains unproved for unbounded labels.

The parity-chain edge `b=1` is genuinely separate.  The analogous pair of
one-sided edge bounds fails, as does the stronger proposal to inject every
`(a,1)` source only into the two global multiplicity spaces of labels `a-1`
and `a+1`.  The latter already fails for plus word `[2,2,2,2,2,2,3]` and
target `(4,1)`, with deficit `410`.  Hence invariant complementary boundary
components are essential at the edge; a complete induction still needs an
edge carrier switch in addition to `(OSD)`.

The C++ program `character_ring_iter/search_su2_defect_cone.cpp` searches
this stronger inequality using exact integers.  The step-one diagonal
monotonicity `E(a,b)>=E(a-1,b-1)` is false: plus labels `[1,3,3,3]` give a
difference `-1` at `(3,3)`, and even after disjoint-support reduction plus
labels `[2,2,2,2]` give a difference `-5` at `(4,4)`.  The parity-respecting
two-layer inequality passed all 1,716 sorted plus multisets with labels at
most six and at most seven factors after imposing the required disjoint
support condition.  This computation motivates `(TLS)` but does not prove
it.

The still stronger claim that multiplication by a disjoint plus label has a
nonnegative matrix in the individual parity-ray basis is also false.  The
exact verifier `character_ring_iter/verify_su2_tls_operator.cpp` finds

```text
plus label 1, target (2,2), input ray (1,4): coefficient -1.
```

Thus the compensation in `(TLS)` groups several rays.

A second tempting localization, to a fixed assignment of the plus factors
between the two tensor coordinates, is false as well.  For a bipartition
`P=A disjoint-union B`, write `m_A(r)` and `m_B(r)` for the corresponding
ordinary tensor-product multiplicities.  Pairing an assignment with its
coordinate reversal gives the proposed local summand

```text
[m_A(a+b)+m_A(a+b-2)]m_B(0)
  + m_A(a-2)m_B(b-2) - m_A(a)m_B(b)
  + (A exchanged with B).
```

The exact C++ diagnostic
`character_ring_iter/search_su2_tls_bipartition.cpp` finds the first negative
value at plus word `[1,1,1,1,1,3]`, split as
`[1,1,1,1] disjoint-union [1,3]`, and target `(a,b)=(2,2)`: the paired
contribution is `-2`.  The complete `(TLS)` packet for that word and target
is nevertheless `48`.  Its nonzero oriented assignment profile is

```text
first side []                total   4
first side [1,1]             total -10
first side [1,1,1,1]         total -15
first side [1,3]             total   5
first side [1,1,1,3]         total  40
first side [1,1,1,1,1,3]     total  24
```

Therefore the required switching map must move plus factors between the two
coordinates; neither fixed-bipartition positivity nor a direct sum of its
coordinate-reversal pairs can prove `(TLS)`.

## The exact `USp(4)` packet form of two-layer switching

Let `Psi_(r,s)`, `r>=s>=0`, denote the irreducible `USp(4)` character of
highest partition `(r,s)`.  Its restriction to the common maximal torus can
be written entirely in `SU(2)` characters:

```text
Psi_(r,s)(x,y)
  = [chi_(r+1)(x) chi_s(y) - chi_s(x) chi_(r+1)(y)]/(x-y).
```

In particular `D_a=Psi_(a-1,0)`.  The normalized `USp(4)` Weyl measure is

```text
dnu(x,y) = (1/2)(x-y)^2 dmu(x)dmu(y),
```

because `integral (x-y)^2 dmu dmu=2`.  Consequently

```text
E_P(a,b) = integral_USp(4) D_a D_b F_P dnu,
F_P(x,y) = product_(p in P) [chi_p(x)+chi_p(y)].
```

The following identities isolate exactly which `USp(4)` coefficients are
tested by `(TLS)`.

**Proposition 8 (two-layer packet identity).** Assume `a>=b>=1`, use
`D_0=0`, and set

```text
R_(a,b) = D_a D_b - D_(a-2)D_(b-2).
```

Then `R_(a,b)` is character-positive and has the multiplicity-free expansion

```text
R_(a,b)
 = sum_(i=0)^(b-1) Psi_(a+b-2-i,i)
   + sum_(i=1)^(b-1) Psi_(a+b-3-i,i-1).          (P8.1)
```

Moreover every plus generator has the three-term virtual-character formula

```text
chi_p(x)+chi_p(y)
  = Psi_(p,0) - Psi_(p-1,1) + Psi_(p-2,0),       (P8.2)
```

where a character with an invalid index is zero.

**Proof.** First, for `m>=n>=0`, the rank-two symplectic row-product rule is

```text
Psi_(m,0) Psi_(n,0)
  = sum_(0<=j<=i<=n) Psi_(m+n-i-j,i-j).          (P8.3)
```

This follows directly from the displayed torus formula, without importing a
tensor-product table.  Multiply `(P8.3)` by `x-y`.  On the left use

```text
(x-y)Psi_(m,0) = chi_(m+1)(x)-chi_(m+1)(y),
Psi_(n,0) = sum_(t=0)^n chi_t(x)chi_(n-t)(y).
```

In the first product apply the `SU(2)` Clebsch--Gordan rule to
`chi_(m+1)chi_t`.  Its terms are indexed by `(t,q)`, `0<=t<=n` and
`0<=q<=t`.  The bijection

```text
i=n-t+q,  j=q
```

maps them exactly to `0<=j<=i<=n`; the terms with `x` and `y` exchanged give
the second half of each alternating numerator.  This proves `(P8.3)`.

Apply `(P8.3)` with `(m,n)=(a-1,b-1)`.  The terms in the product with
`(m,n)=(a-3,b-3)` are exactly those obtained by replacing `(i,j)` with
`(i+2,j+2)`.  Subtraction therefore leaves precisely the layers `j=0` and
`j=1`, which are the two sums in `(P8.1)`.

For `(P8.2)`, multiply its right-hand side by `x-y` and use the torus formula.
The result is

```text
chi_(p+1)(x)-chi_(p+1)(y)
-chi_p(x)chi_1(y)+chi_1(x)chi_p(y)
+chi_(p-1)(x)-chi_(p-1)(y).
```

The recurrence `chi_(p+1)+chi_(p-1)=chi_1 chi_p` reduces this expression to
`(x-y)[chi_p(x)+chi_p(y)]`, proving `(P8.2)`.  QED.

If

```text
F_P = sum_(r>=s>=0) f_P(r,s) Psi_(r,s)
```

is its virtual `USp(4)` expansion, orthonormality and Proposition 8 give the
exact packet formula

```text
Q_P(a,b)
 = sum_(i=0)^(b-1) f_P(a+b-2-i,i)
   + sum_(i=1)^(b-1) f_P(a+b-3-i,i-1).           (P8.4)
```

Thus the remaining `(TLS)` statement is not positivity of every virtual
character coefficient in `(P8.2)`: that is false already for one generator.
It is positivity of the two adjacent coefficient packets `(P8.4)` whenever
`a` and `b` are absent from the plus-label support.

**Corollary 9 (a character-positive sufficient cone).** If the virtual
`USp(4)` expansion of `F_P` is character-positive, then `(TLS)` holds for
every target.  In particular, it holds for a plus word consisting of any
number of fundamental labels `1`.

**Proof.** Every coefficient in each packet `(P8.4)` is nonnegative under the
hypothesis.  For label `1`, formula `(P8.2)` reduces to

```text
chi_1(x)+chi_1(y) = Psi_(1,0).
```

Therefore a word of `n` fundamental plus labels has
`F_P=Psi_(1,0)^n`, an actual tensor-power character, and hence has a
character-positive irreducible expansion.  QED.

Combined with Proposition 7, this proves the complete ordinary two-minus
sector whenever all plus factors are fundamental.  This is an infinite
family, but it does not cover a general plus word because `(P8.2)` has a
negative middle term for every label at least two.

There is also a reproducible smoothing phenomenon beyond Corollary 9.  The
exact C++ program `character_ring_iter/search_su2_sp4_smoothing.cpp` expands

```text
Psi_(1,0)^n [chi_p(x)+chi_p(y)]
```

in irreducible `USp(4)` characters by the coefficient identity

```text
f_(r,s)
 = c_(r,s) + c_(r+2,s) - c_(r+1,s-1) - c_(r+1,s+1).
```

For every `p=2,...,12`, the first character-positive power is exactly
`n=p^2-3`, and positivity persists through `n=150`.  The observed thresholds
are

```text
p                 2   3   4   5   6   7   8   9  10  11  12
first n           1   6  13  22  33  46  61  78  97 118 141
```

Once a power is character-positive, persistence for all later powers is
automatic because multiplication by the genuine character `Psi_(1,0)`
preserves character positivity.  The formula `p^2-3` is currently an exact
finite observation, not a proved general threshold.  It supplies a possible
oscillating-tableau route to a larger sufficient cone, but does not establish
`(TLS)` for unsmoothed general words.

That route has the following exact non-computational reduction.  Let
`K_n(lambda,mu)` denote the number of length-`n` nearest-neighbor walks from
`lambda+(2,1)` to `mu+(2,1)` in the strict type-`C_2` chamber

```text
x_1 > x_2 > 0,
```

with step set `{+e_1,-e_1,+e_2,-e_2}`.  Equivalently, `K_n(lambda,mu)` is
the multiplicity of `Psi_mu` in `Psi_(1,0)^n Psi_lambda`.

**Proposition 10 (reflected-walk form of one-label smoothing).** The
coefficient of `Psi_lambda` in

```text
Psi_(1,0)^n [chi_p(x)+chi_p(y)]
```

is exactly

```text
K_n((p,0),lambda)
  - K_n((p-1,1),lambda)
  + K_n((p-2,0),lambda).                         (P10.1)
```

If `U_n(d_1,d_2)` is the number of unrestricted coordinate-step walks of
length `n` and displacement `(d_1,d_2)`, put

```text
eta_0=(p+2,1), eta_1=(p+1,2), eta_2=(p,1),
xi=lambda+(2,1).
```

Then `(P10.1)` is the fully explicit alternating sum

```text
sum_(w in W(C_2)) det(w)
  [U_n(xi-w eta_0) - U_n(xi-w eta_1) + U_n(xi-w eta_2)].   (P10.2)
```

Here `W(C_2)` is the eight-element group of signed coordinate permutations,
and

```text
U_n(d_1,d_2)
 = sum_(j=0)^n binom(n,j)
     binom(j,(j+d_1)/2)
     binom(n-j,(n-j+d_2)/2),                    (P10.3)
```

with a binomial coefficient interpreted as zero when its lower argument is
not an integer in range.

**Proof.** Formula `(P8.2)` gives the three initial characters in `(P10.1)`.
The minuscule tensor rule for the defining `USp(4)` representation says that
multiplication by `Psi_(1,0)` moves a dominant partition by one of
`+e_1,-e_1,+e_2,-e_2` whenever the result remains dominant.  After the Weyl
shift by `(2,1)`, this is precisely a walk in the displayed strict chamber,
proving `(P10.1)`.

Because the steps have coordinate length one, a path can leave the chamber
only by first hitting one of its two walls.  Reflect the initial segment at
the first wall hit.  The standard sign-reversing reflection pairs all paths
from the signed Weyl images of the start except the paths that stay strictly
inside.  This proves the type-`C_2` reflection formula used in `(P10.2)`.
Finally, choose the `j` coordinate-one steps, then choose their positive
steps and the positive coordinate-two steps.  This gives `(P10.3)`.  QED.

Thus the conjectural threshold `n>=p^2-3` is now a definite eight-image
binomial inequality.  Proving that inequality would enlarge Corollary 9 to
every word with one label `p` and at least `p^2-3` fundamental labels.  It
would still be a sufficient subcone, not the full all-character theorem.

## The exact wreath-product path target

There is a genuine representation behind the nonnegative coefficient matrix
`C_P`, even though its `USp(4)` index is virtual.  Let

```text
K = (SU(2) x SU(2)) semidirect C_2,
```

where the nontrivial element exchanges the two factors, and let `E_p` be the
representation induced from `V_p tensor 1`.  Its restriction to the identity
component has character

```text
chi_p(x)+chi_p(y).
```

Consequently `F_P` is the restriction of the genuine `K`-representation
`tensor_(p in P) E_p`.

An ordered fusion path for a plus word `P=(p_1,...,p_N)` starts at `(0,0)`.
At step `j`, it chooses one coordinate and replaces its current label `u` by
any label

```text
r in u star p_j = {|u-p_j|, |u-p_j|+2, ..., u+p_j},
```

leaving the other coordinate fixed.  The coefficient `c_P(a,b)` is exactly
the number of these paths ending at `(a,b)`.

**Proposition 11 (path form of two-layer switching).** For `a,b>=2`, `(TLS)`
is precisely the path-count inequality

```text
#Path_P(a,b)
 <= #Path_P(a-2,b-2)
      + #Path_P(a+b,0) + #Path_P(a+b-2,0).       (P11.1)
```

Moreover, if `a` and `b` are absent from the label support of `P`, then in
every path counted on the left, the final update of each coordinate starts
from a nonzero label.

**Proof.** The plus-factor update is

```text
C_(pP) = N_p C_P + C_P N_p,
```

which is exactly the stated choice of a coordinate followed by one
Clebsch--Gordan transition.  Induction from `C_empty=e_0 e_0^T` proves the
path interpretation.  Substitution in `(TLS)` gives `(P11.1)`.  Finally, a
Clebsch--Gordan transition out of label zero has the unique output `p_j`.
If the last update producing `a` started at zero, then `a=p_j`, contrary to
the support assumption; the argument for `b` is identical.  QED.

Each transition has a canonical refined-block encoding.  If

```text
r = u+p-2q,  0<=q<=min(u,p),
```

encode it by `q` down-steps followed by `p-q` up-steps.  This block never
goes below zero and uniquely records the Clebsch--Gordan output.  Lowering
the final output by two amounts to replacing `q` by `q+1`; it is possible
unless the final block is saturated at `q=min(u,p)`.  Thus a prospective
injection proving `(P11.1)` is already canonical on paths whose two final
blocks are unsaturated.  The only unresolved paths are those with at least
one saturated final block.  They must be switched into the two boundary
families ending at `(a+b,0)` and `(a+b-2,0)` or paired with another
assignment.  The fixed-bipartition counterexample above proves that this
last switching necessarily changes which coordinate receives some plus
factors.

The still stronger claim that multiplying a packet by one disjoint plus
generator is `USp(4)`-character-positive is false.  The exact C++ diagnostic
`character_ring_iter/search_su2_tls_packet_product.cpp` finds the first
interior failure at plus label `3`, target `(2,2)`, and character `(3,2)`,
with coefficient `-1`.  Hence the path switch cannot be replaced by ordinary
character positivity of each one-step packet product.

The refined blocks isolate a smaller residual than `(P11.1)`.  Order the plus
labels nondecreasingly.  Call the last transition made on a coordinate
*saturated* if its output is the smallest Clebsch--Gordan output, equivalently
if its block parameter is `q=min(u,p)`.  Let `U_P(a,b)` count paths to `(a,b)`
for which the last transitions on both coordinates are unsaturated.

**Proposition 12 (canonical lower-layer injection).** If `a,b>=2`, then

```text
U_P(a,b) <= c_P(a-2,b-2).                        (P12.1)
```

Consequently `(TLS)` follows from the residual saturation estimate

```text
c_P(a,b)-U_P(a,b)
  <= c_P(a+b,0)+c_P(a+b-2,0).                   (SAT)
```

**Proof.** In a path counted by `U_P(a,b)`, increase the contraction parameter
`q` by one in the final block on each coordinate.  Unsaturation says that
`q+1<=min(u,p)`, so both modified transitions are valid.  Their outputs are
`a-2` and `b-2`, and no later transition uses either output.  The assignments,
all earlier channels, and the two modified blocks determine the source path,
so this is an injection and proves `(P12.1)`.  Adding `(SAT)` gives

```text
c_P(a,b)
 <= c_P(a-2,b-2)+c_P(a+b,0)+c_P(a+b-2,0),
```

which is `(TLS)`.  QED.

The exact flagged-path program
`character_ring_iter/search_su2_tls_saturation.cpp` tested `(SAT)` for every
nondecreasing plus multiset with labels at most ten and at most nine factors:
all 92,378 words passed for every disjoint target.  This is strictly
stronger evidence than the earlier direct `(TLS)` check, because it verifies
the proposed decomposition into the explicit lower-layer injection and the
saturated residual.

The order is essential.  On the same path set in decreasing label order,
`(SAT)` fails for plus word `[3,3,3,2,1,1,1]` and target `(5,5)`: the saturated
residual is `126`, while the two boundary families have size `121`.  Thus the
remaining proof target is the following precise lemma:

```text
For nondecreasing plus labels and targets absent from their support,
inject every path having a saturated final block into the union of the
two boundary path families ending at (a+b,0) and (a+b-2,0).
```

Processing the labels from small to large is substantive: the final factor
assigned to either coordinate is then the largest factor assigned to that
coordinate.  Any successful saturation switch must use this extremality;
an order-free final-block argument is ruled out by the displayed descending
counterexample.

This terminology is not merely analogous to symmetric-chain theory.  Let
`C_u={0,...,u}` and `C_p={0,...,p}` be ranked chains.  If `u>=p`, decompose
their rectangle into the chains indexed by `q=0,...,p`:

```text
(0,q),(1,q),...,(u-q,q),(u-q,q+1),...,(u-q,p).
```

The chain starts in rank `q`, ends in rank `u+p-q`, and therefore has length
label `u+p-2q`.  The construction with the coordinates exchanged covers
`p>u`.  These chains are disjoint and exhaust the rectangle, so this is a
symmetric-chain decomposition of `C_u x C_p`.  Its chain labels are exactly
the Clebsch--Gordan outputs, and its innermost chain `q=min(u,p)` is exactly
the saturated output used above.  Iterating the rectangle construction gives
the fusion-path basis.  Thus `(SAT)` is equivalently a matching theorem for
innermost chains in the two-colour product decomposition.

The general connection between the de Bruijn--Tengbergen--Kruyswijk
construction, products of chains, and `sl_2` Jordan chains is developed by
Srinivasan, *Symmetric chains, Gelfand--Tsetlin chains, and the Terwilliger
algebra of the binary Hamming scheme*
(<https://arxiv.org/abs/1001.0280>).  That theory validates the chain model,
but its stated product theorem does not supply the two-colour boundary
matching `(SAT)`; the latter remains the new lemma required here.

One tempting localization of this matching is false.  Partition a path by
the two indices at which its coordinates are updated for the last time.  Since
the final word index is one of them, write the pair as `{t,N}`.  One could ask
for the saturated residual in each fixed `t` class to be bounded by the
boundary paths with the same pair of last-update indices.  The exact C++
diagnostic `character_ring_iter/search_su2_tls_last_update.cpp` disproves this:
for the ascending plus word `[1,1,1,3,3]`, `t=3`, and target `(5,2)`, the
fixed-class residual is `6` while the corresponding boundary count is `0`.
Therefore a successful injection for `(SAT)` must be allowed to change the
last-update index as well as the coordinate assignment.  This is compatible
with the earlier fixed-bipartition obstruction and narrows the target to a
global carrier-switching map across the ordered word.

There is nevertheless a stronger cumulative statement which survives that
obstruction.  Write the ordered word as

```text
P=(p_1,...,p_N),       p_1<=...<=p_N.
```

Partition a residual path by the earlier of its two last-update indices.  Let
`R_i(a,b)` be the saturated residual in the class `{i,N}`.  Define `B_i(r)`
analogously for boundary paths ending at `(r,0)`, and let `B_empty(r)` count
the boundary paths in which the zero coordinate is never updated.  Exact
enumeration supports the prefix Hall inequalities

```text
sum_(i<=t) R_i(a,b)
 <= B_empty(a+b)+B_empty(a+b-2)
      +sum_(i<=t)[B_i(a+b)+B_i(a+b-2)]           (PH_t)
```

for every `t<N`.  The final prefix is exactly `(SAT)`.  Unlike the false
pointwise comparison, `(PH_t)` has an explicit cut interpretation.  Let
`C_t` be the two-colour coefficient matrix after the prefix
`(p_1,...,p_t)`, and put

```text
K_t=N_(p_(t+1)) ... N_(p_N).
```

Then last-update decomposition gives the exact identity

```text
B_empty(r)+sum_(i<=t) B_i(r)
  = (K_t C_t+C_t K_t)_(r,0).                    (CUT_t)
```

Indeed, after the cut all remaining factors must update one coordinate: the
two terms in `(CUT_t)` record which coordinate it is.  The never-updated
case is included in the zero column of `C_t`.  Moreover

```text
K_t=sum_r m_(p_(t+1),...,p_N)(r) N_r
```

has nonnegative integral coefficients.  Thus `(PH_t)` says that the
saturated part of the same monochromatic-suffix transfer is dominated by
the two outer boundary channels.  This is the precise algebraic form of a
leftward carrier switch; it also explains why boundary supply created at an
earlier last-update time can pay a later deficit.

The refined strict C++ diagnostic
`character_ring_iter/search_su2_tls_last_update.cpp --prefix` checked
`(PH_t)` for all 1,716 nondecreasing words with labels at most six and at
most seven factors, every disjoint target, and every prefix.  Its minimum
prefix slack was zero only in a vacuous zero-count case.  Among prefixes with
positive residual, the minimum slack was five, at plus word `[1,1,3,3]`,
cut `t=3`, and target `(4,2)`, where the cumulative residual and boundary
counts were respectively four and nine.  This is stronger evidence for a
Hall matching, but it is not a proof of `(PH_t)`.

A larger exact run has now checked all 12,870 nondecreasing words with labels
at most eight and at most eight factors.  Every prefix inequality passed.
This run used the same exact integer C++ implementation and strengthens only
the finite evidence, not the unbounded Hall lemma.

The same prefix formulation includes the parity-chain edge `b=1`.  Under the
disjoint-support hypothesis every plus label is different from `1`.  If a
transition with label `p` ends at `1`, Clebsch--Gordan parity forces its input
to be `p-1` or `p+1`; hence the transition is automatically saturated.  Thus
the residual at `(a,1)` is the entire source family, while the two boundary
channels in `(PH_t)` are exactly `a+1` and `a-1`.  The final prefix inequality
is therefore the edge estimate

```text
c_P(a,1) <= c_P(a+1,0)+c_P(a-1,0).
```

The `--prefix-edge` mode checks only this extension.  Through labels at most
six and at most seven factors it passed all 1,716 ordered words and every
prefix.  The smallest positive-residual slack was `13`, for plus word
`[2,2,3,3]`, cut `t=3`, and target `(5,1)`, with residual `4` and cumulative
boundary supply `17`.  Hence one proof of the extended prefix-Hall lemma
would simultaneously provide the missing edge carrier switch and the
interior saturation estimate.

There is an exact exterior-square formulation which makes the edge carrier
literal.  Let `V` be the free module with basis `e_0,e_1,...`, put

```text
A_p=N_p tensor I+I tensor N_p,
D=N_1 tensor I-I tensor N_1.
```

All fusion matrices commute, so `D A_p=A_p D`.  The plus coefficient vector
is

```text
C_P=product_(p in P) A_p (e_0 tensor e_0),
```

and

```text
D(e_0 tensor e_0)=e_1 tensor e_0-e_0 tensor e_1.
```

Taking the coefficient of `e_a tensor e_0` therefore gives the exact
identity

```text
d_P(a-1)+d_P(a+1)-c_P(a,1)
 = coefficient_(e_a wedge e_0)
     [product_(p in P) dGamma(N_p)](e_1 wedge e_0),       (EDGE-EXT)
```

where `dGamma(N)=N tensor I+I tensor N` restricted to the exterior square.
At finite level the same identity holds with the missing upper-wall neighbor
omitted.  Thus the edge problem, ordinary and finite-level, is precisely
positivity of one matrix coefficient of an ordered product of additive
compounds.  The initial `e_1 wedge e_0` is the proposed `B_1` carrier.

Equivalently, introduce independent variables and set

```text
K(z)=product_i (I+z_i N_(p_i)).
```

Then the right side of `(EDGE-EXT)` is the square-free coefficient

```text
[z_1...z_n] det K(z)[rows (a,0), columns (1,0)].          (EDGE-MINOR)
```

Indeed, exterior powers turn the product defining `K` into the product of
the exterior powers of its factors, and selecting every variable once
selects `dGamma(N_(p_i))` from each factor.  This identifies the missing
edge theorem as a coefficientwise order-two total-positivity statement,
but only for this square-free coefficient and these boundary indices.

The prefix-Hall strengthening has the same form.  For

```text
K_t=N_(p_(t+1)) ... N_(p_N),
```

its edge slack is the matrix coefficient from `e_1 wedge e_0` to
`e_a wedge e_0` of

```text
(K_t tensor I+I tensor K_t)
  product_(i<=t) dGamma(N_(p_i)).                         (EDGE-CUT)
```

Thus `(PH_t)` asks for positivity after the remaining suffix has been forced
to travel as a single monochromatic carrier.  This is a sharper target for a
reflection or carrier-switching proof than the raw two-colour count.

The strict exact C++ verifier
`character_ring_iter/search_su2_level_edge_exterior.cpp` implements
`(EDGE-CUT)` directly in every finite fusion ring.  It passed all 48,619
nondecreasing plus words with labels in `2,...,k`, at most eight factors, and
levels `2<=k<=10`, checking every cut and every target absent from the word.
This is the first cut-by-cut evidence that the same carrier cone survives the
affine upper wall; it is finite evidence, not the required uniform proof.

The ordinary edge admits a still sharper two-packet decomposition.  Write

```text
X_p=e_(p+1) wedge e_0,
Y_p=e_(p-1) wedge e_0-e_p wedge e_1.             (PACKET)
```

For `p>=2`, the first exterior update is exactly

```text
dGamma(N_p)(e_1 wedge e_0)=X_p+Y_p.             (P12B.1)
```

If `Q=(q_1,...,q_m)` and `T_Q=product_i dGamma(N_(q_i))`, then the outer
packet has the non-circular interpretation

```text
coefficient_(e_a wedge e_0) T_Q X_p=E_Q(a,p+1). (P12B.2)
```

Indeed, under the `USp(4)` identification the target row is `D_a` and
`X_p=D_(p+1)`, so `(P12B.2)` is the defining defect integral.  It is therefore
nonnegative under induction on the number of plus factors.  The only new
edge target is the following ordered packet lemma:

```text
If 2<=p<=q_1<=...<=q_m, then
coefficient_(e_a wedge e_0) T_Q Y_p>=0
for every a.                                      (OPD)
```

Combining `(P12B.1)`, `(P12B.2)`, and `(OPD)` proves the ordinary `b=1` edge
by induction: choose the smallest plus label as `p`, use the earlier
two-minus theorem on the outer packet, and use `(OPD)` on the other packet.
No disjoint-target qualification is needed inside `(OPD)`.

There is also an exact defect-only form of the new packet coefficient:

```text
coefficient_(e_a wedge e_0) T_Q Y_p
 = sum_(r in a star p) E_Q(r,1)+E_Q(a,p-1)
     -E_Q(a-1,p)-E_Q(a+1,p).                   (P12B.3)
```

Thus `(OPD)` is precisely the compensation missing from the false
four-term Monge inequality above.  The optimized strict C++ diagnostic
`character_ring_iter/search_su2_edge_packet_domination.cpp` checked both
packets first through labels `2,...,8` and eight factors, and then through
labels `2,...,10` and nine factors.  The larger run enumerated 48,620
multisets (the empty bookkeeping word and 48,619 nonempty words) and passed
every target, including targets belonging to the plus support.  The size
order is substantive: singling out the largest packet instead fails for word
`[2,2,2,2,3,3,4]`, target `5`, with ordered packet coefficient `-17`.

A fixed-colour proof of `(OPD)` is impossible.  Discrete integration by
parts in the internal split `j=1,...,p-1` reduces the slack for a fixed
suffix bipartition with fusion kernels `K,L` to

```text
K(a,p-1)L(0,0)+K(a,1)L(0,p)
 -K(a,0)L(0,p-1)-K(a,p)L(0,1).                (P12B.4)
```

The exact diagnostic `search_su2_edge_packet_split.cpp` finds value `-5`
for word `[2,2,2,2,3]`, the split assigning the whole suffix to `L`, and
target `2`.  Even pairing a split with its complement fails: word
`[2,2,2,3,3]`, suffix mask `4`, and target `3` give `-4`.  Hence the proof
of `(OPD)` must switch factor assignments globally; local recoupling within
each fixed colour is insufficient.

Parity does not reduce that global switch to a single odd carrier.  The
strict diagnostic `character_ring_iter/search_su2_edge_packet_toggle.cpp`
pairs assignments differing only by the side receiving the first odd suffix
factor.  It fails for

```text
plus=[5,5,6,6,6,6,6],  target=5,  paired slack=-424.
```

Even after summing over every assignment of all even suffix factors, while
fixing only the other odd-factor colours, the first failure is

```text
plus=[5,5,5,6,6,6,6],  target=6,  grouped slack=-176.
```

This obstruction persists for the smallest packet `p=2`.  The restricted
mode of the same strict C++ diagnostic finds

```text
plus=[2,7,8,8,8,8,8,8], target=2, paired slack=-4916,
plus=[2,7,7,8,8,8,8,8], target=3, even-summed slack=-2400.
```

The second line already sums every assignment of the even factors.

Thus a mixed-parity proof has to redistribute several odd factors and the
even factors globally; neither a one-carrier toggle nor an even-factor
average is a sufficient local block.

There is a useful dual form of `(OPD)`.  Keep `X_r=Psi_(r,0)` and define

```text
Z_(a,p)=X_(a-1)[Psi_(p-1,1)-X_(p-2)]
          + indicator_(a<p) X_(p-a-1).                 (P12C.1)
```

**Proposition 12C (multiplicity-free dual packet).** `Z_(a,p)` is a genuine
multiplicity-free `USp(4)` character for every `a>=1`, `p>=2`, and

```text
X_(a-1)Y_p
 = indicator_(a<p) X_(p-a-1)-Z_(a,p).                  (P12C.2)
```

Its support is explicit.  It is the disjoint union of the following two
families of partitions `(r,s)`:

```text
A_(a,p): r=a+s-p-1,  1<=s<=p-1,  a>=p+1;

B_(a,p): r=a+1-p+2u-s,
         1<=u<=p-1,  1<=s<=u<=r,
         and either r>u or u=p-1.                       (P12C.2a)
```

In particular,

```text
Z_(a,p)=sum_((r,s) in A_(a,p) union B_(a,p)) Psi_(r,s). (P12C.2b)
```

**Proof.** The second identity follows immediately from
`Y_p=X_(p-2)-Psi_(p-1,1)`.  For positivity use the rank-two symplectic Pieri
rule in its elementary witness form: the multiplicity of `Psi_mu` in
`X_m Psi_lambda` is the number of partitions `nu` such that both
`lambda/nu` and `mu/nu` are horizontal strips and

```text
|lambda/nu|+|mu/nu|=m.                                  (P12C.3)
```

Here compare `lambda_0=(p-2,0)` with `lambda_1=(p-1,1)` and put `m=a-1`.
Every Pieri witness for `lambda_0` is `nu=(t,0)`, `0<=t<=p-2`.  Map it to a
witness for `lambda_1` as follows.  If `mu_1>=t+1`, use
`nu'=(t+1,0)`.  Horizontality of `mu/(t,0)` gives `mu_2<=t`, so the new
skew shape is horizontal, and increasing `|nu|` by one on both sides keeps
`(P12C.3)` unchanged.  If `mu_1=t` and `mu_2>=1`, instead use
`nu'=(t,1)`; the same size identity and horizontality are immediate.

The only witness not mapped occurs when `mu=(t,0)`.  Then `(P12C.3)` forces

```text
t=p-a-1,
```

so this is the single negative row constituent `X_(p-a-1)`, and it exists
exactly when `a<p`.  The map is injective.  Conversely, an unmapped witness
for `lambda_1` is of one of two types.  A witness `(u,0)` is unmapped exactly
when `mu_2=u`; its size equation gives family `A_(a,p)`.  A witness `(u,1)`
is unmapped when `mu_1>u`, and also in the upper-boundary case `u=p-1`, where
the prospective `lambda_0` witness does not exist.  Its size equation gives
family `B_(a,p)`.  The upper-boundary case is essential; for example it
produces `Psi_(3,2)` when `(p,a)=(4,2)`.

The size equation fixes `u` in either type.  Both types cannot occur for the
same `mu`: their fixed values differ by one, while horizontality of the
second type requires `u>=mu_2`.  There are likewise no repetitions within a
type.  Thus `(P12C.2a)`--`(P12C.2b)` follow and every positive excess
multiplicity is exactly one.  Adding the row correction in `(P12C.1)`
cancels the unique negative constituent and leaves a multiplicity-free
genuine character.  QED.

If

```text
F_Q=sum_(r>=s) f_Q(r,s) Psi_(r,s),
```

then `(P12C.2)` turns ordered packet domination into the single cumulative
coefficient inequality

```text
<Z_(a,p),F_Q>
 <= indicator_(a<p) f_Q(p-a-1,0),                       (P12C.4)
```

for `p<=min(Q)`.  In particular, for `a>=p` the sum of the coefficients of
`F_Q` over a definite multiplicity-free dual packet must be nonpositive.
The strict C++ program `inspect_su2_edge_packet_dual.cpp --verify 30`
independently checked the complete support formula `(P12C.2a)`, the zero/one
expansion, and the unique row correction in all 870 pairs `2<=p<=30`,
`1<=a<=30`.

The cumulative inequality also has a formulation involving only genuine
`SU(2) x SU(2)` representations.  Put `K=SU(2) x SU(2)`, write `V_r` for
the `SU(2)` irreducible of label `r`, and define

```text
M_p = sum_(j=1)^(p-1) V_j tensor V_(p-j),
W_(a,p;Q) = Res_K(X_(a-1)) tensor M_p
              tensor product_(q in Q)
                [(V_q tensor 1) direct-sum (1 tensor V_q)].       (P12D.1)
```

Every representation in `(P12D.1)` is genuine, and `W_(a,p;Q)` is invariant
under exchange of the two factors.  Let `m_(r,s)` denote the multiplicity of
`V_r tensor V_s` in this representation.

**Proposition 12D (low-`K`-type form of ordered packet domination).**  For
every `a>=1`, `p>=2`, and suffix `Q`, one has the exact identity

```text
coefficient_(e_a wedge e_0) T_Q Y_p
   = m_(1,1)-m_(0,0)-m_(2,0).                            (P12D.2)
```

Consequently `(OPD)` is equivalent to the genuine multiplicity inequality

```text
m_(1,1) >= m_(0,0)+m_(2,0)                              (P12D.3)
```

for the modules `(P12D.1)` with `p<=min(Q)`.  Thus it is enough to inject
the direct sum of the `(0,0)` and `(2,0)` multiplicity spaces into the
`(1,1)` multiplicity space.

**Proof.**  The row-character branching formula and `(P8.2)` give

```text
Res_K(X_r)=sum_(i=0)^r V_i tensor V_(r-i),
Res_K(Y_p)=-M_p.                                         (P12D.4)
```

By the normalized measure formula preceding Proposition 8, the left side
of `(P12D.2)` is

```text
-(1/2) integral_K W_(a,p;Q)(x,y) (x-y)^2 dmu(x)dmu(y).
```

Here `x=chi_1(x)` and `y=chi_1(y)`, so the Clebsch--Gordan rule gives

```text
(x-y)^2
 = 2 + chi_2(x)+chi_2(y)-2 chi_1(x)chi_1(y).             (P12D.5)
```

Taking the trivial-`K` multiplicity in `(P12D.5)` times `W`, and using
exchange symmetry `m_(0,2)=m_(2,0)`, makes the integral equal to
`2[m_(0,0)+m_(2,0)-m_(1,1)]`.  This proves `(P12D.2)` and `(P12D.3)`.  QED.

This form exposes the exact hierarchy needed for induction.  If `W` is any
exchange-symmetric genuine `K`-module with multiplicities `m_(r,s)` and

```text
W' = W tensor [(V_q tensor 1) direct-sum (1 tensor V_q)],
```

then, for `q>=2`, Clebsch--Gordan gives

```text
m'_(1,1) = 2[m_(q-1,1)+m_(q+1,1)],
m'_(0,0) = 2m_(q,0),
m'_(2,0) = m_(q-2,0)+m_(q,0)+m_(q+2,0)+m_(q,2).         (P12D.6)
```

Thus the next OPD index is exactly

```text
2m_(q-1,1)+2m_(q+1,1)
 -m_(q-2,0)-3m_(q,0)-m_(q+2,0)-m_(q,2).                 (P12D.7)
```

Formula `(P12D.7)` is a concrete ordered low-type matching target.  A proof
still has to exploit that `W` starts with the first two factors in
`(P12D.1)` and that every previously adjoined suffix label lies between `p`
and `q`; arbitrary symmetric genuine `K`-modules do not satisfy it.
The independent strict C++ verifier
`character_ring_iter/inspect_su2_edge_packet_kindex.cpp` constructs the
genuine module `(P12D.1)` directly in the `K` representation ring and
compares `(P12D.2)` with the exterior-square packet coefficient.  It passed
all 1,694 cases through labels eight, two suffix factors, and every supported
target in that domain.  The same run independently checked the closed
two-suffix formula used in Proposition 12E.

The low index cannot be embedded in a uniform diagonal-concavity cone.  For
`(p,a;Q)=(2,4;[2,3])`, the desired `r=1` slack is `3`, but the formally
analogous next slack is

```text
m_(2,2)-m_(1,1)-m_(3,1)=-19.
```

Hence `(P12D.3)` is a genuinely bottom-of-chamber inequality; an induction
requiring the same diagonal inequality at every height is false.

The first nontrivial part of `(OPD)` can in fact be proved uniformly.

**Proposition 12E (two-suffix ordered packet domination).**  Ordered packet
domination holds when the suffix `Q` has at most two factors.

**Proof.**  The empty suffix is immediate from the definition of `Y_p`.
For one factor `q>=p`, only the first summand of `Y_p` can reach the boundary,
apart from the orientation reversal `e_0 wedge e_1=-e_1 wedge e_0` when
`q=p`.  Hence

```text
coefficient_(e_a wedge e_0) dGamma(N_q)Y_p
 = indicator_(a in q star (p-1))
    +indicator_(q=p) indicator_(a=1) >= 0.              (P12E.1)
```

Now take `p<=q<=r`, put `W=dGamma(N_q)Y_p`, and use symmetry of the fusion
matrices.  The final boundary update is

```text
[dGamma(N_r)W]_(a,0)
 = sum_(t in a star r) W_(t,0)+W_(a,r).                 (P12E.2)
```

The first term in `(P12E.2)` is the cardinality of

```text
(a star r) intersect (q star (p-1)).                    (P12E.3)
```

plus the nonnegative correction
`indicator_(q=p) indicator_(1 in a star r)`.

A direct expansion of the four wedges in `dGamma(N_q)Y_p` gives

```text
W_(a,r)
 = indicator_(r=q) indicator_(a=p-1)
   +indicator_(a=1) indicator_(r in q star p)
   -indicator_(a=p) indicator_(r in q star 1)
   +indicator_(r=p) indicator_(a in q star 1).          (P12E.4)
```

All terms in `(P12E.4)` are nonnegative except the third.  Under
`p<=q<=r`, that term is nonzero only when `a=p` and `r=q+1`.  In this case

```text
q star (p-1) is contained in p star (q+1)
```

and the smaller fusion interval has exactly `p` elements.  Thus `(P12E.3)`
contributes `p`, the unique negative term contributes `-1`, and the total is
at least `p-1>0`.  This proves the two-factor case.  QED.

The next suffix length can also be closed uniformly.  The proof below is
computer-assisted only in one quantifier-free Presburger lemma; all
representation-theoretic reductions and the exceptional family are
explicit.

**Proposition 12E1 (three-suffix ordered packet domination).**  Ordered
packet domination `(OPD)` holds when the suffix has three factors.  More
precisely, if

```text
2<=p<=q<=r<=s,
```

then

```text
coefficient_(e_a wedge e_0)
  dGamma(N_s)dGamma(N_r)dGamma(N_q)Y_p >= 0             (P12E1.1)
```

for every `a>=1`.

**Proof.**  For a submultiset `A` of `L={q,r,s}`, write `m_A(t)` for the
multiplicity of `V_t` in the tensor product of the labels in `A`.  Expanding
which of the three suffix factors updates each exterior carrier gives the
exact eight-term identity

```text
C_a=sum_(B subset L) [
       m_B(0)m_((L minus B) union {p-1})(a)
      -m_B(p-1)m_(L minus B)(a)
      -m_B(1)m_((L minus B) union {p})(a)
      +m_B(p)m_((L minus B) union {1})(a)].             (P12E1.2)
```

The `B=empty` term is

```text
M=m_(p-1,q,r,s)(a).
```

All positive terms other than `M` will be denoted by `P`.  Since every
suffix label is at least `p`, the complete negative load in `(P12E1.2)` is

```text
N=sum_({u,v,w}=L) [
      indicator_(p-1 in u star v) indicator_(a=w)
     +indicator_(1 in u star v) indicator_(a in w star p)]
   +indicator_(a=p)m_(q,r,s)(1).                        (P12E1.3)
```

Thus `C_a=M+P-N` with `P>=0`.

Use the fusion tree which first couples `p-1` with `q`, then couples the
result with `r`, and finally with `s`.  It gives

```text
M=#{(x,y): x in (p-1) star q,
             y in x star r,
             a in y star s}.                           (P12E1.4)
```

There is a fixed finite witness pool inside this path set.  Put

```text
x_i=q-p+1+2i,                          0<=i<=5,
l_i=max(|x_i-r|,|a-s|),
y_(i,j)=l_i+2j,                        0<=j<=5,
```

and retain `(x_i,y_(i,j))` exactly when

```text
x_i<=q+p-1,
x_i+r and a+s have the same parity,
y_(i,j)<=min(x_i+r,a+s).                           (P12E1.5)
```

Every retained pair is a distinct path in `(P12E1.4)`.  The following finite
Presburger statement is exact:

```text
#retained pairs >= N,                                 (P12E1.6)
```

except when

```text
p is even,       q=r=p,       s=p+1,       a=p.        (P12E1.7)
```

The strict C++/Z3 verifier
`character_ring_iter/verify_su2_opd_three_suffix_z3.cpp` encodes
`2<=p<=q<=r<=s`, the fusion-interval indicators in `(P12E1.3)`, and all 36
candidate predicates in `(P12E1.5)`.  It asks for a nonexceptional integer
solution with fewer candidates than `N`; Z3 4.8.12 proves the negation
unsatisfiable.  No label cutoff occurs in this query.  Therefore `M>=N` and
hence `C_a>=0` outside `(P12E1.7)`.

It remains to check the exceptional family.  Write `p=2n`.  In the pairing
tree `(p-1,p)|(p,p+1)`, the main multiplicity counts

```text
0<=i<=2n-1,   0<=j<=2n,
|i-j|<=n,     i+j>=n-1.
```

Consequently

```text
M=5n(n+1)/2.                                          (P12E1.8)
```

Indeed, the band `|i-j|<=n` contains `3n^2+2n` pairs, and the excluded
lower triangle contains `n(n-1)/2`.  In `(P12E1.2)` the six negative unit
terms consist of two `p-1` pair channels, two adjacent pair channels, and
the two channels contributing to `m_(p,p,p+1)(1)`.  The positive remainder
also equals six: the two singleton `p` terms contribute two each, while the
equal pair `(p,p)` contributes two more.  Hence `P=N=6` and `C_a=M>0`.
This proves `(P12E1.1)`.  QED.

The independent strict C++ program
`character_ring_iter/analyze_su2_opd_three_suffix.cpp` constructs the
exterior updates directly, checks `(P12E1.2)`, enumerates the 36-channel
witness pool, and verifies the exceptional formula.  Through label 15 it
passed 85,680 ordered parameter/target cases.  Its most negative interior
packet is `-6`; the boundary reservoir is then at least `11`.

The subset-cardinality grading gives a useful exact reformulation of the
general packet, but two natural filtrations of it are false.  Put

```text
M_p=sum_(j=1)^(p-1) V_j tensor V_(p-j),
W_Q(z)=M_p product_(i=1)^n
       [(V_(q_i) tensor 1)+z(1 tensor V_(q_i))].        (P12E2.1)
```

If `W_j` is the coefficient of `z^j`, and `(W_j)_b` is its multiplicity
space for the second `SU(2)` type `V_b`, then the cardinality-`j` contribution
to the boundary coefficient is the first-factor virtual character

```text
f_j=(W_j)_1-chi_1(W_j)_0.                              (P12E2.2)
```

Indeed, under the alternating-character map

```text
e_r wedge e_s |-> chi_r(x)chi_s(y)-chi_s(x)chi_r(y),
```

one has

```text
Y_p |-> (y-x)M_p.
```

Each suffix update multiplies this character by
`chi_(q_i)(x)+chi_(q_i)(y)`.  Extracting the coefficient of
`chi_a(x)chi_0(y)` gives `(P12E2.2)`, coefficient by coefficient in `z`.
In particular,

```text
coefficient_(e_a wedge e_0) T_QY_p
 = coefficient_(chi_a) sum_(j=0)^n f_j.                (P12E2.3)
```

Thus one could try to prove the apparently stronger one-sided prefix claim

```text
sum_(j=0)^h f_j belongs to the nonnegative character cone.              (PREF_h)
```

Exact enumeration initially supported this through six suffix factors and
labels at most six, but `(PREF_h)` is false.  The first counterexample found
in the enlarged run is

```text
p=2, Q=[2,2,3,3,4,4,4,4], a=1,
(f_0,...,f_8)=(1456,1044,22,-620,-1348,-564,598,612,1456).
```

The prefix through `h=5` is `-10`, while the full ordered-packet coefficient
is `2656`.  Hence the full theorem genuinely uses assignments beyond a
one-sided cardinality cut.

Pairing low-cardinality layers with their coordinate-reversed high-cardinality
layers does not repair the filtration.  The exchange-symmetric shell claim

```text
sum_(j<=h or j>=n-h) f_j belongs to the nonnegative character cone       (SHELL_h)
```

fails for

```text
p=2, Q=[3,3,3,3,3,3], a=3,
(f_0,...,f_6)=(210,0,600,-520,390,-240,0).
```

The shell through depth one is `210+0-240+0=-30`, although the full
coefficient is `440`.  Therefore a proof of general `(OPD)` must mix several
interior cardinality layers; neither Boolean prefix Hall matching nor an
outside-in complement matching can be the missing uniform injection.

The same calculation identifies the correct unfiltered index.  Since
`W_Q=W_Q(1)` is exchange-symmetric,

```text
coefficient_(e_a wedge e_0) T_QY_p
 = <chi_a tensor chi_1,W_Q>
    -<chi_(a-1) tensor 1,W_Q>
    -<chi_(a+1) tensor 1,W_Q>.                          (P12E2.4)
```

Equivalently this is the row coefficient of the `USp(4)` virtual character

```text
[Psi_(p-2,0)-Psi_(p-1,1)]
 product_(q in Q)[Psi_(q,0)-Psi_(q-1,1)+Psi_(q-2,0)]   (P12E2.5)
```

against `Psi_(a-1,0)`.  Formula `(P12E2.5)` is an unfiltered
Dirac-induction/Pieri target.  Standard symplectic Pieri positivity does not
apply directly because every bracketed suffix factor is virtual.  The
ordered hypothesis `p<=min(Q)` must therefore be used in a cancellation or
straightening theorem for the row coefficient itself.

The strict arbitrary-precision C++ tools
`character_ring_iter/search_su2_opd_subset_prefix.cpp` and
`character_ring_iter/inspect_su2_opd_subset_prefix.cpp` reproduce both
counterexamples.  The searcher's `prefix` and `shell` modes check the two
failed filtrations separately.

There is nevertheless an unbounded dominant-largest-label extension.  It
uses the full exterior carrier rather than either false cardinality cut.

**Proposition 12E3 (dominant suffix transfer).**  Let

```text
R=(q_1,...,q_m),       p<=q_1<=...<=q_m,
T=q_1+...+q_m,
W=T_RY_p.
```

Assume every boundary coefficient of `W` is nonnegative.  If `q>=q_m`
(with the evident convention for empty `R`) and

```text
q>=p+T-3,                                                   (P12E3.1)
```

then every boundary coefficient of `dGamma(N_q)W` is nonnegative.

**Proof.**  Write `b_t` for the coefficient of `e_t wedge e_0` in `W`, and
extend the wedge coefficients antisymmetrically in their two indices.  The
last update gives the exact identity

```text
[dGamma(N_q)W]_(a,0)
 =sum_(t in a star q)b_t+W_(a,q).                          (P12E3.2)
```

The first term is nonnegative by hypothesis.  Every wedge in `W` has index
sum at most

```text
M=p+1+T
```

and index-sum parity `p-1+T`.  If `q>p+T`, this support bound makes
`W_(a,q)=0` for every `a>=1`.  If `q=p+T`, only `a=1` can occur.  The unique
top-degree path starts in `-e_p wedge e_1`, sends every suffix factor through
the first coordinate, and takes every outer Clebsch--Gordan branch.  Hence

```text
W_(q,1)=-1,                 W_(1,q)=1.                    (P12E3.3)
```

Finally let `q=p+T-1`.  Support and parity leave only `a=2`.  A top-degree
path ending at the unordered pair `{q,2}` exists only when `p=2`: every
suffix update then travels through the second coordinate and takes its outer
branch.  Its original orientation is `-e_2 wedge e_q`, so

```text
W_(2,q)=-indicator_(p=2).                                (P12E3.4)
```

In the exceptional case `p=2`, the boundary coefficient `b_q` is exactly
one.  This can be seen directly from the four terms in the subset formula:
the empty second-colour subset and the top branch from `e_1 wedge e_0`
give one; every negative term would require a nonempty second-colour subset
of total label at most one, impossible because all suffix labels are at
least two; the same support bound excludes any additional positive term.
Since `q` belongs to `2 star q`, `(P12E3.2)` contains this unit `b_q`, which
cancels `(P12E3.4)`.  All remaining terms in this case are nonnegative.

It remains in fact to handle one more value allowed by `(P12E3.1)`.  Let
`q=p+T-2`, put `m=|R|`, and let `n_2` be the number of labels equal to two
in `R`.  Now `M=q+3`, so support and parity leave only `a=1,3`.  Counting
the unique total-defect-two location in the two tensor coordinates gives

```text
W_(1,q)=m+indicator_(p=2)(1+n_2),
W_(3,q)=indicator_(p=2)n_2-indicator_(p=3).              (P12E3.5)
```

For completeness, the `m` terms in the first line are the choices of the
single next-to-outer Clebsch--Gordan branch in the product with initial
label `p`.  When `p=2`, assigning all suffix factors to the other coordinate
gives the additional unit, and assigning one label-two factor there gives
the additional `n_2` choices.  At target three, a label-two factor gives
the `n_2` positive paths for `p=2`; the only negative path occurs for
`p=3`, when all suffix factors travel through the second coordinate.

In that last case the boundary coefficient `b_(q+1)` is exactly one.  The
top path from `e_(p-1) wedge e_0` gives the unit.  In the subset formula any
negative contribution would require a nonempty second-colour subset of
total label at most one, while every such label is at least three; all other
positive terms are excluded by the same support bound.  Since `q+1` belongs
to `3 star q`, this unit cancels the final `-1` in `(P12E3.5)`.  This proves
this case.

For the last value `q=p+T-3`, support and parity leave only `a=2,4`.  Let
`n_3` count the labels equal to three in `R`.  The same outer-branch count
now gives

```text
W_(2,q)=
  n_3-2n_2-m,                 p=2;
  1+n_3,                      p=3;
  0,                          p>=4,

W_(4,q)=
  n_3-n_2,                    p=2;
  n_3,                        p=3;
  -1,                         p=4;
  0,                          p>=5.                       (P12E3.6)
```

These formulas again just locate either the unique total-defect-two branch
or a top branch whose coordinate receives a subset of total label at most
four.  The only deficit of unbounded size is the first line with `p=2`.
In that case direct extraction from the four subset terms gives the three
outer boundary coefficients

```text
b_(q+2)=1,
b_q=m+n_2,
b_(q-2)=m(m+1)/2-1+n_2(n_2-1)+n_2(m-1)-n_2 n_3.         (P12E3.7)
```

All three belong to `2 star q`.  Adding them to the first line of
`(P12E3.6)` gives

```text
m(m+1)/2+n_2(n_2+m-n_3-3)+n_3.                          (P12E3.8)
```

This is nonnegative.  If `n_2=0` it is immediate.  If `n_2=1`, write
`m=1+n_3+r`, `r>=0`; the possibly negative second term is `r-1`, and the
whole expression is still nonnegative, with equality only at
`(m,n_2,n_3)=(1,1,0)`.  If `n_2>=2`, then
`n_2+m-n_3-3=2n_2+r-3>=1`.

At target four with `p=2`, just `b_(q+2)+b_q` and the first line for
`W_(4,q)` already sum to `1+m+n_3`.  For `p=3` both carrier entries are
nonnegative.  For `p=4`, the unique `-1` at target four is cancelled by the
top boundary coefficient `b_(q+2)=1`, and `q+2` belongs to `4 star q`.
This exhausts `(P12E3.6)` and proves the proposition.  QED.

Consequently a minimal counterexample to `(OPD)`, under induction on suffix
length, must satisfy the balanced-support condition

```text
q_m<=p+q_1+...+q_(m-1)-4.                                (P12E3.9)
```

Combining Proposition 12E1 with repeated dominant suffix transfer proves
`(OPD)` for arbitrary-length ordered sequences satisfying

```text
q_j>=p+q_1+...+q_(j-1)-3       for every j>=4.            (P12E3.10)
```

Thus the arbitrary-length theorem is now complete on a superincreasing
mixed-parity cone as well as on the all-even cone.  The strict exact C++
verifier `character_ring_iter/verify_su2_opd_dominant_suffix.cpp` constructs
the exterior state, checks `(P12E3.3)`--`(P12E3.7)` and every exceptional
reservoir, and then checks every extended boundary.  Through prefix labels
seven and seven prefix factors it passed 3,002 prefixes, 17,988 dominant
extensions, and 1,120,764 boundary targets.

There is also a complete parity subcase.  Let

```text
F_Q=sum_(r,s>=0) c_Q(r,s) chi_r(x)chi_s(y)
```

be the genuine `K`-character expansion, so `c_Q(r,s)` counts the ordered
two-colour fusion paths from `(0,0)` to `(r,s)`.  Multiplying the two source
wedges in `Y_p` gives the exact boundary formula

```text
coefficient_(e_a wedge e_0) T_QY_p
 = sum_(r in a star (p-1)) c_Q(r,0)-c_Q(a,p-1)
    -sum_(r in a star p) c_Q(r,1)
    +sum_(r in a star 1) c_Q(r,p).                       (P12F.1)
```

**Proposition 12F (all-even ordered packet domination).**  If `p` and every
label in `Q` are even, then `(OPD)` holds for every target `a`.

**Proof.**  Starting from `(0,0)`, assigning only even labels to either
coordinate can produce only even endpoint labels.  Hence

```text
c_Q(a,p-1)=0,                 c_Q(r,1)=0
```

because `p-1` and `1` are odd.  The two remaining terms in `(P12F.1)` are
sums of genuine path counts and are nonnegative.  QED.

The parity proof has a stronger invariant-cone form which isolates the next
mixed sector.  For `u>v` of opposite parity define

```text
B_(u,v)= e_u wedge e_v,          u odd,
B_(u,v)=-e_u wedge e_v,          u even.                (P12F.2)
```

**Proposition 12F1 (even-prefix signed wedge cone).**  If `p` is even and
every label in `R` is even, then

```text
T_RY_p=sum_(u>v, u-v odd) h_(u,v) B_(u,v),
h_(u,v)>=0.                                             (P12F.3)
```

**Proof.**  Initially

```text
Y_p=B_(p-1,0)+B_(p,1).
```

It is enough to check one even update.  In a term where `N_(2r)` updates
the larger index `u`, every output has the parity of `u`.  If it remains
larger than `v`, the orientation and the sign in `(P12F.2)` are unchanged.
If it crosses below `v`, wedge orientation contributes a minus sign, while
the new larger index has the opposite parity and the basis sign contributes
the same minus sign.  Their product is positive.  The identical argument
applies when the smaller index is updated and crosses above the larger one.
Thus `dGamma(N_(2r))` has a nonnegative integral matrix in the basis
`(P12F.2)`, proving `(P12F.3)` by induction.  QED.

Conceptually, this cone is the image of an ordinary character cone.  Apply
central translation in the second `SU(2)` factor,

```text
J(chi_r tensor chi_s)=(-1)^s chi_r tensor chi_s.
```

For opposite-parity `u,v`, the alternating character of `B_(u,v)` is sent
to the genuine symmetric character
`chi_u tensor chi_v+chi_v tensor chi_u`.  Every even plus generator is fixed
by `J`.  This gives a second proof of `(P12F.3)` and explains why the cone is
special to an all-even prefix.

The twist also connects this sector exactly to the fundamental-difference
normal form below.  Put

```text
S_r=chi_r tensor 1+1 tensor chi_r,
H_p=J(Y_p)=S_(p-1)+chi_p tensor chi_1+chi_1 tensor chi_p.
```

The Clebsch--Gordan recurrence gives the positive-seed identity

```text
H_p=S_p S_1-S_(p+1).                                  (P12F.3a)
```

For `p=2` it factors further as

```text
H_2=(chi_1 tensor chi_1)S_1.                           (P12F.3b)
```

If `H=H_p product_(r in R)S_r=sum_b H_b tensor chi_b` and `q` is odd, then
central translation changes its plus update into a difference update:

```text
J(T_qT_RY_p)
 =(chi_q tensor 1-1 tensor chi_q)H
 =Delta D_q H.                                         (P12F.3c)
```

Here the last equality is Proposition 23.  Extracting the second-factor
trivial type in `(P12F.3c)` gives the representation-valued identity

```text
(id tensor tau)J(T_qT_RY_p)=chi_q H_0-H_q.              (P12F.3d)
```

Thus `(EOH)` is equivalently the one-minus fundamental-transfer assertion

```text
chi_q H_0-H_q belongs to the nonnegative character cone                 (EOH')
```

for the structured seed `(P12F.3a)` followed by even `S_r` generators.
This is much smaller than the false assertion for arbitrary symmetric
positive `H`: it retains exactly the generator grouping required by the
counterexample after Proposition 23.

There is an exact two-reservoir decomposition of `(EOH')`.  Write

```text
C=product_(r in R)S_r=sum_b C_b tensor chi_b.
```

Parity and `(P12F.3a)` give

```text
H_0=chi_(p-1)C_0+chi_1 C_p,
H_q=sum_(b in q star (p-1))C_b
       +chi_p(C_(q-1)+C_(q+1)).                         (P12F.3e)
```

Consequently

```text
chi_qH_0-H_q
 =sum_(b in q star (p-1))(chi_b C_0-C_b)
   +chi_p[(chi_(q-1)+chi_(q+1))C_p
             -C_(q-1)-C_(q+1)].                        (P12F.3f)
```

The first line is a sum of ordinary boundary-defect characters for the
all-even plus word `R`; in an induction proving the two-minus theorem it is
therefore supplied by the outer hypothesis.  The second line carries the
new internal-column correction.

Neither reservoir can be separated further.  The internal-column packet is
already negative for

```text
p=2, q=3, R=[4], target=2:
first packet=2, internal packet=-1, total=1.             (P12F.3g)
```

Nor does the first packet alone dominate the negative half of the internal
packet: for `p=2`, `q=3`, `R=[2,2,2]`, and target two, that proposed
reservoir has value `-10`.  Hence the positive term
`chi_p(chi_(q-1)+chi_(q+1))C_p` in `(P12F.3f)` is essential.  The strict
arbitrary-precision C++ diagnostic
`character_ring_iter/analyze_su2_eoh_decomposition.cpp` evaluates
`(P12F.3f)` and searches these subpackets.  Through even labels ten, eight
even factors, and odd `q<=13`, it checked 2,001 words and 11,010 `(p,q;R)`
cases: the complete packet was always nonnegative, the first defect packet
was nonnegative throughout, and both proposed separations failed as above.

The same target has a closed one-step recurrence.  For any symmetric
character `K=sum_b K_b tensor chi_b`, define

```text
F_b(K)=chi_b K_0-K_b.
```

Then a direct use of the fusion rule gives

```text
F_q(K S_r)
 =chi_r F_q(K)+sum_(s in q star r)F_s(K)-chi_q F_r(K).   (P12F.3h)
```

Indeed, `(K S_r)_0=chi_rK_0+K_r` and
`(K S_r)_q=chi_rK_q+sum_(s in q star r)K_s`; substituting and using
`chi_qchi_r=sum_(s in q star r)chi_s` proves the identity.  Thus the exact
inductive closure statement for `(EOH')` is the structured fusion-concavity
inequality

```text
chi_r F_q(K)+sum_(s in q star r)F_s(K)>=chi_q F_r(K),    (EFC)
```

where `q` is odd, `r` is an even next label at least `p`, and
`K=H_p product S_(r_i)`.  The negative final term explains both failed
reservoir separations.  Proving `(EFC)` by a Temperley--Lieb differential
would make the one-odd cone invariant under every further even update.

This turns the case of exactly one odd suffix label into a purely unsigned
matching problem.  Let `q` be odd, let `W=T_RY_p`, and extend `h` symmetrically
to unordered pairs.  The last-update identity and `(P12F.3)` give, for even
`a`,

```text
[dGamma(N_q)W]_(a,0)
 =sum_(t in a star q) h_(t,0)-h_(a,q),                  (P12F.4)
```

while for odd `a` both sides vanish by parity.  Consequently OPD with even
`p`, an all-even prefix, and one odd suffix is equivalent to the unsigned
endpoint inequality

```text
h_(a,q)<=sum_(t in a star q)h_(t,0).                    (EOH)
```

Unlike the false fixed-colour and odd-toggle statements, `(EOH)` has no
signed summands: `h` counts paths in the positive wedge graph generated by
`(P12F.2)`.  A proof still needs an injection from paths ending at the
unordered pair `{a,q}` to boundary paths whose odd endpoint lies in
`a star q`.

The strict exact C++ verifier
`character_ring_iter/verify_su2_opd_even_prefix_cone.cpp` checks the signed
cone, the transfer identity `(P12F.4)`, and `(EOH)`.  Through even labels
ten, eight even prefix factors, and odd labels at most thirteen it passed
2,001 prefixes, 612,841 signed-cone entries, 11,010 odd extensions, and
592,566 boundary targets.  This is exact finite evidence for `(EOH)`, not
its uniform proof.

The unsigned injection still cannot be made atom by atom.  Even under the
correct label order, the strict diagnostic
`character_ring_iter/search_su2_even_wedge_boundary.cpp --ordered-atoms`
finds

```text
p=2, q=5, initial atom=e_6 wedge e_2,
even suffix=[2,2,2,6,6,6], target=2, boundary coefficient=-3.
```

The complete four-atom packet in this case has coefficient `795`.  Thus a
proof of `(EOH)` must mix the initial atoms as well as the subsequent colour
assignments.  The same diagnostic also shows why dropping the order is
unsafe: the unrestricted initial wedge `e_6 wedge e_4` followed by
`[2,2,2,2,2,6]` has boundary coefficient `-25`.

For the full edge packet, rather than the lower packet `Y_p` alone, central
translation gives a stronger wall-compatible reduction.  It is uniform in
the ordinary ring and every finite fusion ring.

**Proposition 12F2 (parity-twist edge normal form).**  Let `E` be the
multiset of even plus labels and `O` the multiset of odd plus labels.  In the
ordinary ring or `SU(2)_k`, put `B_1=e_1 wedge e_0`.  Then

```text
J(T_(E union O)B_1)
 =Delta^|O| S_1 product_(o in O)D_o product_(e in E)S_e. (P12F.5)
```

In particular, if `O={q}` and `C=product_(e in E)S_e=sum_b C_b tensor e_b`,
then

```text
(id tensor tau)J(T_qT_EB_1)
 =sum_(s in q star_k 1)(e_s C_0-C_s).                   (P12F.6)
```

Here `star_k` means finite fusion, and in the ordinary ring it is the usual
Clebsch--Gordan interval.  Consequently the full edge defect satisfies the
exact reduction

```text
E_(E union {q})(a,1)
 =sum_(s in q star_k 1) E_E(a,s).                       (P12F.7)
```

**Proof.**  One has `J(B_1)=S_1`.  An even additive update is fixed by `J`,
whereas an odd update changes sign on the second tensor factor:

```text
J dGamma(N_e) J^(-1)=e_e tensor 1+1 tensor e_e=S_e,
J dGamma(N_o) J^(-1)=e_o tensor 1-1 tensor e_o=Delta D_o.
```

All factors commute, and Proposition 23 supplies the last equality in the
second line, proving `(P12F.5)` ordinarily and at finite level.  For one odd
label, `C` has only even-even types.  Hence

```text
(S_1C)_0=e_1C_0,
(S_1C)_q=sum_(s in q star_k 1)C_s.
```

Also `e_qe_1=sum_(s in q star_k 1)e_s`, with the upper neighbor simply
absent at the affine wall.  Extracting the second-factor trivial type from
`(e_q tensor 1-1 tensor e_q)S_1C` proves `(P12F.6)`.  The coefficient of
`e_a` in `e_sC_0-C_s` is exactly `E_E(a,s)`, giving `(P12F.7)`.  QED.

Thus, in the induction on the number of plus factors, every edge word with
exactly one odd plus label is already closed: the right side of `(P12F.7)`
contains only defects for the strictly shorter all-even word `E`.  This does
not require separate positivity of `Y_p`.

That distinction is essential at finite level.  The strict exact verifier
`character_ring_iter/verify_su2_level_even_prefix_cone.cpp` finds

```text
k=5, p=2, even prefix=[2,2,2,2,4,4,4], q=5, a=2:
Y_p packet=-5, outer packet=28, full edge=23.             (P12F.8)
```

The lower-packet EOH therefore does not survive the wall, while the combined
identity `(P12F.7)` does.  The independent strict C++ verifier
`character_ring_iter/verify_su2_level_one_odd_edge_reduction.cpp` checked
247,246 instances of `(P12F.7)` through level twelve and seven even factors.
The packet verifier checked 294,700 boundary targets in the same level range:
the lower packet failed as displayed, but every combined edge packet passed.

For two odd plus labels the same calculation leaves one genuinely coupled
fusion-concavity inequality.

**Proposition 12F3 (two-odd coupled edge reduction).**  Let `q,r` be odd,
let `E` be an even plus word, put

```text
C=product_(e in E)S_e,
F_b(C)=e_b C_0-C_b.
```

Then, ordinarily and in `SU(2)_k`, the full edge character is

```text
(id tensor tau)J(T_qT_rT_EB_1)
 =e_q sum_(v in r star_k 1)F_v(C)
   +e_r sum_(u in q star_k 1)F_u(C)
   -e_1 sum_(s in q star_k r)F_s(C).                   (P12F.9)
```

**Proof.**  Put `K=S_1C`.  Since `C` has only even-even types,

```text
K_0=e_1C_0,
K_q=sum_(u in q star_k 1)C_u,
K_r=sum_(v in r star_k 1)C_v,
K_s=e_1C_s                    for even s.
```

Extracting the second-factor trivial type from
`(e_q tensor 1-1 tensor e_q)(e_r tensor 1-1 tensor e_r)K` gives

```text
e_qe_rK_0-e_qK_r-e_rK_q+sum_(s in q star_k r)K_s.
```

Substitute the preceding four identities and then replace
`C_b=e_bC_0-F_b(C)`.  Every multiple of `C_0` cancels by the fusion identity;
the remaining terms are exactly `(P12F.9)`.  QED.

**Corollary 12F3A (bottom odd terminals are removed by support
reduction).**  In an induction on the number of signed factors, the
two-odd edge target `(P12F.9)` may be restricted, ordinarily and in every
`SU(2)_k`, to

```text
q,r>=3.                                                (P12F.9a)
```

Equivalently, after the odd-label sign toggle `(P5A.1)`, neither of the two
minus terminals can equal the mandatory plus label `f=1` in a
support-disjoint minimal case.

**Proof.**  If, say, `q=1` in the original form `(P12F.9)`, the signed word
contains both the plus factor `s_1` and the edge minus factor `d_1`.  At
level `k>=2`, the fusion rule `1 star_k 1={0,2}` gives

```text
s_1 d_1=d_0+d_2=d_2,
```

whereas at level `k=1` it gives `s_1d_1=d_0=0`.  Thus the former case is a
signed word with one fewer factor and the latter vanishes.  The induction
hypothesis proves nonnegativity.  The argument for `r=1` is identical.
This is also the specialization `i=1` of Proposition 5, and it is
independent of the number and labels of the even background factors.  QED.

The last line of `(P12F.9)` cannot be paid by either positive family alone.
The strict exact verifier
`character_ring_iter/verify_su2_level_two_odd_edge_reduction.cpp` finds the
minimal separate failures

```text
k=3, E=empty, q=1, r=3, a=3: second family minus sink=-1;
k=5, E=[2,2], q=r=3, a=5: first family minus sink=-1.
```

Through level twelve and seven even factors it checked 819,454 instances of
`(P12F.9)`.  The identity held exactly and the complete coupled expression
was always nonnegative.  This remains finite evidence: the uniform proof now
requires one map from the sink
`e_1 sum_(s in q star_k r)F_s(C)` into the direct sum of both positive
families.

The earlier two-layer switching cone by itself does not supply that map.
Suppose abstractly that a symmetric even-index defect matrix has zero
boundary and nonnegative diagonal increments

```text
Q(x,y)=E(x,y)-E(x-2,y-2).
```

Substituting `E(x,y)=sum_t Q(x-2t,y-2t)` into the coefficient form of
`(P12F.9)` produces a negative ray even in the parity-toggle-minimal regime:

```text
a=5, q=r=3, a+1=q+r,
coefficient of the diagonal-increment ray based at (6,6) = -1.          (P12F.10)
```

The strict combinatorial C++ verifier
`character_ring_iter/verify_su2_two_odd_tls_cone.cpp` computes every ray
coefficient directly.  Thus symmetry, zero boundary, and `(TLS)` do not
imply `(P12F.9)`; a proof must exploit the smaller defect *orbit* generated
by Proposition 6, or the stronger ordered one-sided data `(OSD)`.

There is a sharper four-terminal form of the same obstruction.  First apply
Proposition 5A to all four odd labels.  At the level of defects this gives
the exact identity

```text
E_(E union {q,r})(a,1)=E_(E union {a,1})(q,r).          (P12F.11)
```

Thus `(P12F.9)` is the two-minus inequality with minus terminals `q,r`,
plus terminals `a,1`, and an all-even background.  In the subset expansion,
every nonzero invariant block contains zero, two, or four of these odd
terminals.  The pairings `{q,r}|{a,1}` have positive sign, while the four
cross-pairings have negative sign.  This is the precise four-terminal
Temperley--Lieb problem; it is not merely an analogy.

A fixed allocation of the even background is not nonnegative.  The strict
finite-fusion diagnostic
`character_ring_iter/search_su2_four_terminal_packets.cpp` finds, already
in the stable range,

```text
k=20, E=[2,2], q=r=3, a=5, one label 2 on each side:
zero/four-terminal part=0, positive pairings=0,
crossed pairings=4, signed packet=-4.
```

For the four allocations of the two indexed even factors the complete
profile is

```text
17, -4, -4, 17,       total=26.                        (P12F.12)
```

Hence any switching or Plucker map must move even factors between the two
tensor coordinates.  Uncrossing inside a fixed background split cannot
prove `(P12F.9)`.

The exterior state gives an exact source--curl normal form.  For opposite
parities put

```text
M_(u,v)=e_u tensor e_v+e_v tensor e_u.
```

Then, ordinarily and in every `SU(2)_k`, one has

```text
J(T_q T_r B_1)
 =S_1 sum_(s in q star_k r) S_s
   -sum_(u in q star_k 1) M_(u,r)
   -sum_(v in r star_k 1) M_(v,q).                    (P12F.13)
```

**Proof.**  On the exterior square, commutativity and the fusion rule give

```text
T_q T_r B_1
 =sum_(s in q star_k r)T_s B_1
   +sum_(u in q star_k 1)e_u wedge e_r
   +sum_(v in r star_k 1)e_v wedge e_q.
```

Every `s` is even, so central translation commutes with `T_s` and sends
`T_sB_1` to `S_sS_1`.  Every `u,v` is even whereas `q,r` are odd; central
translation sends each displayed cross wedge to minus the corresponding
`M`.  This proves `(P12F.13)`.  The argument uses only fusion associativity
and parity, so the affine-wall truncations are already included.  QED.

After multiplication by `C_E=product_(e in E)S_e`, `(P12F.9)` is therefore
equivalent to the source--curl Hall inequality

```text
[e_a tensor e_0] C_E S_1 sum_(s in q star_k r)S_s
 >=[e_a tensor e_0] C_E[
      sum_(u in q star_k 1)M_(u,r)
       +sum_(v in r star_k 1)M_(v,q)].                (SCH)
```

All even updates have nonnegative matrices on the parity-twisted symmetric
basis, but the seed in `(P12F.13)` itself is not in that positive cone.
The strict full-state diagnostic
`character_ring_iter/search_su2_two_odd_twisted_cone.cpp` verifies the
source--curl identity and finds

```text
q=r=3, E=[2]: coefficient of M_(3,2)=-3.              (P12F.14)
```

Thus neither fixed-allocation positivity nor eventual ordinary character
positivity is the missing invariant.  A proof of `(SCH)` must retain a
cumulative dominance relation between the boundary reservoir and the two
adjacent curl bands.  This is now the most economical Temperley--Lieb/Hall
target for the two-odd sector.

The allocation compensation is genuinely nonlocal.  Pairing the two
allocations of a distinguished largest even factor first works in the
profile `(P12F.12)`, but fails in general.  Even summing the entire block of
maximal equal labels is insufficient.  The `--largest-block` mode of the
four-terminal diagnostic finds

```text
k=5, E=[2,2,2,2,2,2,4], q=r=3, a=5:
fixed smaller allocation=[2],
zero/four-terminal part=46, positive pairings=73,
crossed pairings=246, maximal-block orbit=-127.          (P12F.15)
```

For this same case, summing all allocations of a terminal suffix first
becomes nonnegative only when the suffix contains six of the seven even
factors.  Hence no proof which moves only the last factor, only the maximal
label block, or any uniformly bounded terminal block can establish
`(SCH)`.  The correct object must be a cumulative prefix Hall matching of
the kind `(PH_t)`: boundary supply created arbitrarily early in the ordered
word must be allowed to pay a later crossed-pair deficit.

The fundamental carrier can nevertheless be moved exactly inside each
fixed allocation.  Index the even background and, for `T subset E`, put

```text
K_T=product_(e_i in T)N_(e_i),       K_empty=I.
```

Let `h_T(a;q,r)` be the coefficient of `e_a tensor e_0` in `(P12F.13)`
when the first coordinate receives precisely the even factors in `T` and
the second receives those in `T^c`.  Then

```text
h_T(a;q,r)
 =(K_(T^c))_(0,0)(K_T N_q N_r N_1)_(0,a)
   +(K_T)_(1,a)(K_(T^c)N_qN_r)_(0,0)
   -(K_T)_(r,a)(K_(T^c))_(q,1)
   -(K_T)_(q,a)(K_(T^c))_(r,1).                  (P12F.16)
```

**Proof.**  Expand the positive reservoir in `(P12F.13)`.  Its boundary
part gives the first product and its `M_(1,s)` part, summed over
`s in q star_k r`, gives the second.  For the first negative band use

```text
sum_(u in q star_k 1)(K_(T^c))_(u,0)
 =(K_(T^c)N_q)_(0,1)=(K_(T^c))_(q,1),
```

where the middle equality is commutation with the fundamental fusion
matrix.  The other band is identical with `q,r` exchanged.  This proves
the formula.  It is uniform at finite level.  QED.

The full fixed-allocation four-terminal packet is

```text
h_T(a;q,r)+h_(T^c)(a;q,r),                         (P12F.17)
```

and the desired edge coefficient is `sum_T h_T`.  The strict diagnostic
`search_su2_four_terminal_packets.cpp` now verifies `(P12F.16)--(P12F.17)`
independently of its direct sixteen-terminal-assignment calculation.  Thus
the remaining Hall problem has been reduced to four products of entries of
commuting fusion kernels for every allocation.

The proved all-even ordered packet cone does not already contain the
two-odd seed.  For `q=r=3`, exact exterior expansion gives

```text
T_3^2 B_1
 ={1^0:3, 2^1:-1, 3^0:2, 3^2:-2, 4^1:-1,
   4^3:2, 5^0:2, 6^1:-1, 7^0:1}.                 (P12F.18)
```

Here `u^v` denotes `e_u wedge e_v`.  Every degree-seven state of the
controlled forms `T_QX_p`, or `T_QY_p` with even `p<=Q`, has coefficient
at `4^3` at most zero, whereas `(P12F.18)` has coefficient two.  The strict
enumerator `analyze_su2_two_odd_even_packets.cpp` lists the complete finite
set of degree-compatible packets.  Consequently `(P12F.9)` cannot be
deduced by a nonnegative state decomposition into the already-proved
all-even outer and ordered packet cones.

At seven total factors the Plucker multiplication map exposes the first
syzygies in this same family.  For minus labels `[3,3]` and plus labels
`[1,2,2,2,5]`, the 24 negative `3|4` products have rank 22.  A modular
kernel basis consists of two six-term relations, each alternating across
the three repeated label-two positions; the equal-label pair reservoir has
dimension 18.  This agrees with the general kernel target `(7.2)` and shows
that a prefix Hall proof may equivalently be sought as exactness of a
Plucker complex whose first differential resolves these alternating
equal-label relations.

Simultaneous spectral diagonalization does not make the remaining moment
cone pointwise positive.  At level `k`, let `lambda_p(i)` be the eigenvalue
of `N_p` in mode `i`.  The two-odd edge coefficient has the exact spectral
form

```text
sum_(i<j) w_(i,j)(q,r,a)
  product_(e in E)[lambda_e(i)+lambda_e(j)],             (P12F.19)
```

where `w_(i,j)` is the product of the source and target sine minors and the
two odd eigenvalue sums.  For even `e`, the generator coordinate is
unchanged by independently replacing either mode by its central complement,
so one may first sum `w` over the resulting central orbits.

Neither required pointwise sign survives.  The long-double diagnostic
`analyze_su2_two_odd_spectral_orbits.cpp` imposes disjoint support, the
parity-toggle-minimal condition `a+1>=q+r`, and the stable bound
`k>=q+r+a+1`.  Its first case already gives

```text
k=12, q=r=3, a=5:
orbit (1,2), weight 0.0077509557..., label-6 coordinate -1.60387547...;
orbit (2,3), total weight -0.0961509965....               (P12F.20)
```

The margins are far from roundoff, but this is a numerical obstruction,
not an exact sign certificate.  It is sufficient to reject the proposed
proof route: positivity of `(P12F.19)` requires cancellation between
distinct spectral orbits, even in the ordinary stable regime.

The mandatory fundamental plus terminal supplies an explicit Plucker
differential.  This is an exact algebraic statement, not a computational
conjecture.  Let

```text
R_n=k[p_(ij):1<=i<j<=n]/I_(2,n)
```

be the multigraded `SL_2` bracket ring, where `I_(2,n)` is the Plucker
ideal.  Distinguish a vertex `f`, let `R_(n-1)` be the bracket ring on the
other vertices, and let `M_f` be the component of `R_n` having degree one
at `f`.

**Proposition 12F4 (fundamental-column Plucker presentation).**  There is
an exact `R_(n-1)`-module presentation

```text
direct_sum_(i<j<k; i,j,k!=f) R_(n-1) e_(ijk)
  --d--> direct_sum_(i!=f) R_(n-1)e_i
  --phi--> M_f --> 0,                                  (P12F.21)

phi(e_i)=p_(fi),
d(e_(ijk))=p_(jk)e_i-p_(ik)e_j+p_(ij)e_k.              (P12F.22)
```

In a fixed multidegree `p` with `p_f=1`, the coefficient of `e_i` in the
middle term has multidegree `p-e_f-e_i`, while the coefficient of
`e_(ijk)` has multidegree `p-e_f-e_i-e_j-e_k`.

**Proof.**  Every bracket monomial of degree one at `f` contains exactly
one bracket `p_(fi)`, so `phi` is onto.  The Plucker ideal is generated by
its quadratic relations.  Those not involving `f` are already the defining
relations of the coefficient ring `R_(n-1)`.  Every generator involving
`f` has the form

```text
p_(fi)p_(jk)-p_(fj)p_(ik)+p_(fk)p_(ij)=0,
```

and is exactly the image under `phi` of `(P12F.22)`.  Since all relations
linear in the `f`-column are generated over `R_(n-1)` by these quadrics,
`ker(phi)=image(d)`.  The multidegree statement follows term by term.  QED.

This proposition identifies the sought differential in the two-odd edge
sector.  Take `f` to be the mandatory plus label one.  The terminal switch

```text
p_(fq)p_(ra)-p_(fr)p_(qa)+p_(fa)p_(qr)=0              (P12F.23)
```

has exactly the two crossed negative pairings and the positive
`{q,r}|{f,a}` pairing.  Multiplying `(P12F.23)` by arbitrary bracket
polynomials on the even background produces the required prefix-wide
carrier moves; the coefficient is global, so this is consistent with the
unbounded backward depth in `(P12F.15)`.

The ambiguity between such carrier moves is also generated in bounded
degree.  Let `V` be the index space on the vertices other than `f`, grade
`R_(n-1)` by bracket degree, and write

```text
F_0=R_(n-1) tensor V,       F_1=R_(n-1) tensor exterior^3 V.
```

Give `e_i` degree zero and `e_(ijk)` degree one, so the degree-`d` part of
`d:F_1->F_0` is

```text
R_(d-1) tensor exterior^3 V  ->  R_d tensor V.          (P12F.24)
```

**Proposition 12F5 (quadratic generation of carrier syzygies).**  Over a
field of characteristic zero, for every `d>=2`,

```text
ker(d)_d
 = S_(d,d-1,1,1)V direct_sum S_(d-1,d-1,1,1,1)V.     (P12F.25)
```

Moreover the complete `R_(n-1)`-module `ker(d)` is generated by its
degree-two part

```text
ker(d)_2=S_(2,1,1,1)V direct_sum exterior^5 V.         (P12F.26)
```

Thus every relation between prefix-wide Plucker carriers is a bracket
multiple of a polarization of one of two universal moves.  In alternating
index notation they may be represented by

```text
Sigma(x;b,c,d)
 =p_(xb)e_(xcd)-p_(xc)e_(xbd)+p_(xd)e_(xbc),

Tau(x_1,...,x_5)
 =sum_(i<j)(-1)^(i+j+1)p_(x_i x_j)
      e_(x_1 ... omit x_i ... omit x_j ... x_5).       (P12F.27)
```

Here the arguments in `Sigma` are arbitrary vectors, so its polarizations
include the squarefree five-vertex relations, not only relations with one
repeated coordinate index.

**Proof.**  Standard monomial theory gives

```text
R_m = S_(m,m)V.
```

Pieri's rule, with partitions having invalid parts omitted, gives

```text
R_d tensor V
 =S_(d+1,d)V direct_sum S_(d,d,1)V,

R_(d-1) tensor exterior^3 V
 =S_(d,d,1)V
   direct_sum S_(d,d-1,1,1)V
   direct_sum S_(d-1,d-1,1,1,1)V.                    (P12F.28)
```

The map `d` is `GL(V)`-equivariant.  The only common summand in the two
lines is `S_(d,d,1)V`, and `d` is nonzero on it: for example the image of
`p_(12)^(d-1)e_(123)` is nonzero.  Schur's lemma therefore identifies its
image with that summand and proves `(P12F.25)`.

At `d=2` the two kernel summands are exactly `(P12F.26)`.  Direct expansion
of `(P12F.22)` sends `Sigma` to a four-index Plucker relation times `e_x`
and sends `Tau` to a signed sum of the same relations, so both displayed
families lie in the kernel and furnish the two irreducible summands.  For
general `d`, multiply a highest-weight vector in each degree-two summand by
the nonzero highest-weight element `p_(12)^(d-2)` of `R_(d-2)`.  The two
products are nonzero highest-weight vectors of respectively
`S_(d,d-1,1,1)V` and `S_(d-1,d-1,1,1,1)V`.  Since `(P12F.25)` is
multiplicity-free, the `R_(n-1)`-span of `(P12F.26)` contains both complete
summands in every degree.  QED.

The strict multigraded C++ verifier
`character_ring_iter/verify_su2_plucker_second_syzygies.cpp` independently
constructs the full matrix `(P12F.24)` in the noncrossing standard-monomial
basis.  It constructs the pointed four-index moves and all squarefree
five-index polarizations, verifies that their composite with `d` is zero,
and compares their rank with the complete kernel.  The runs

```text
(vertices,maximum coordinate,maximum total degree)
 =(5,3,11): 491 degrees, maximum kernel 22;
 =(6,3,11): 1764 degrees, maximum kernel 70;
 =(7,3,11): 5671 degrees, maximum kernel 146;
 =(8,2,9): 2648 degrees, maximum kernel 161
```

all passed.  These computations audit the formulas and their coordinate
polarizations; the unbounded generation assertion is proved by
`(P12F.28)`, not inferred from the runs.

There are in fact two support kernels, which must not be conflated.  In a
fixed terminal multidegree, strip the unique `f`-edge from each standard
monomial in a negative factorization source.  This gives

```text
iota_-:O_->F_0,                 phi iota_-=mu_-.
```

The restriction of `iota_-` to `ker(mu_-)`, followed by its image, gives
the exact sequence

```text
0 -> ker(iota_-)
  -> ker(mu_-)
  -> image(iota_-) intersect ker(phi) -> 0.            (P12F.29)
```

Surjectivity of the last arrow is just the definition of the intersection.
Since `ker(phi)=image(d)`, the exact kernel count is therefore

```text
kappa_-=kappa_alloc+kappa_car,
kappa_alloc=dim ker(iota_-),
kappa_car=dim(image(iota_-) intersect image(d)).       (P12F.30)
```

The earlier one-stage formulation of `(SCX)` omitted `ker(iota_-)`; the
first seven-factor relations with `a=5` already lie entirely in that
kernel.  A complete support proof must establish two disjoint payments:

```text
kappa_alloc <= P_alloc,       kappa_car <= P_car,
P_alloc+P_car <= P_positive.                          (SCX')
```

The allocation term is the equal-label/transposition problem detected by
the same-sign singlet residues.  The carrier term is the genuine Plucker
problem.  If `L_paid` denotes the carrier preimages supported by the
remaining positive domains, its exact target is

```text
d^(-1)(image(iota_-) intersect image(d))
 subset L_paid + R_(n-1) ker(d)_2.                    (P12F.31)
```

Proposition 12F5 makes `(P12F.31)` local in syzygy type: any two carrier
preimages differ by bracket multiples of `Sigma` and `Tau`.  There are no
additional higher-degree carrier relations to discover.

The strict `--seven-lift` mode of
`verify_su2_plucker_second_syzygies.cpp` constructs `O_-`, `F_0`, `F_1`,
and the global invariant space independently and verifies the commuting
square exactly modulo `1,000,000,007`.  For

```text
[q,r,a,f,E]=[3,3,5,1,2,2,2]
```

it finds `kappa_-=2`, `kappa_alloc=2`, and `kappa_car=0`; both previously
observed six-term relations are allocation relations.  For the doubled
fundamental case

```text
[q,r,a,f,E]=[3,3,1,1,2,2,2]
```

it finds respectively `7`, `2`, and `5`.  The positive `{a,f}` singlet
channel has dimension seven, enough to contain the carrier count, although
an injective unbounded map is not yet proved.

The `--seven-lift-sweep 8` mode checked 300 two-odd edge cases.  Every
commuting square passed.  A nonzero carrier kernel occurred in 23 cases,
all with `a=1`; its largest dimension was five, and in every case it was at
most the dimension of the positive `{a,f}` singlet reservoir.  For `a>1`
the whole negative kernel was already an allocation kernel.  These are
exact modular rank data, not a proof for unbounded labels or factor count.

The earlier `--two-odd-edge` residue test remains relevant to the first
summand: after excluding `a=1`, it passed all 180 cases with labels at most
eight, with maximum allocation-kernel dimension two.  The residual task is
now sharply split: prove the equal-label bound for `kappa_alloc`, and prove
the doubled-fundamental singlet bound for `kappa_car`.  The later
arbitrary-background diagnostic below shows that these seven-factor payments
do not extend one fixed `3|(n-3)` source layer at a time; the global lift must
also use positive supply from other split sizes or the cumulative Hall
formulation.

The doubled-fundamental seven-factor source has a complete elementary
classification.  Write `E={e_1,e_2,e_3}` and, for `s` equal to `q` or `r`,
let `s'` denote the other odd terminal.  Put

```text
U_(s,i)=Inv(V_(s') tensor V_1
             tensor product_(j!=i)V_(e_j)).
```

**Proposition 12F6 (paired orientation of the exceptional source).**  If
`a=f=1`, the negative `3|4` source space is

```text
O_- = direct_sum_(s in {q,r})
        direct_sum_(i: e_i=s-1 or s+1)
          [U_(s,i)^(a) direct_sum U_(s,i)^(f)].        (P12F.32)
```

The superscript records which fundamental lies in the three-element block;
the other fundamental lies in its complement.  In particular the two
orientations have canonically equal dimensions.

**Proof.**  A negative triple contains exactly one of the two minus labels,
say the odd label `s`.  Its other two labels must have odd total parity for
the triple invariant to exist.  Among the plus labels, only `a` and `f` are
odd, while every `e_i` is even.  Hence the triple contains exactly one of
`a,f` and exactly one `e_i`.  The multiplicity-free Clebsch--Gordan rule
gives

```text
Inv(V_s tensor V_1 tensor V_(e_i)) != 0
 iff e_i=s-1 or s+1,
```

and that invariant is then one-dimensional.  The complementary invariant
space is precisely `U_(s,i)`, proving `(P12F.32)`.  QED.

For each summand the two orientations are related by the single bracket
switch

```text
p_(fi)p_(aj)=p_(fa)p_(ij)+p_(fj)p_(ai),               (P12F.33)
```

where `i` and `j` are the respective neighbors of `f` and `a`.  The first
term on the right lies in the positive `{a,f}` singlet channel; the second
is the opposite orientation.  Thus the exceptional carrier bound is now
the following precise paired-switch rank lemma:

```text
dim(image(iota_-) intersect image(d))
 <= dim Inv(V_q tensor V_r tensor V_(e_1)
                         tensor V_(e_2) tensor V_(e_3)).               (PSR)
```

Proposition 12F6 proves that no other source type enters `(PSR)`.  What is
not yet proved is that simultaneous switches from different `(s,i)`
summands have no extra rank loss beyond their common `{a,f}` residual
space.

There is now a concrete generic-rank form of that last assertion.  For each
paired source type `alpha=(s,i,orientation)`, let

```text
rho_alpha:O_- -> Inv(V_q tensor V_r tensor product_j V_(e_j))
```

be zero on all other source types and, on type `alpha`, retain the residual
`p_(ij)` in `(P12F.33)`.  With independent scalars `z_alpha`, put

```text
rho_z=sum_alpha z_alpha rho_alpha.                    (P12F.34)
```

The equal-weight specialization is the global singlet contraction and
vanishes on `ker(mu_-)`.  Generic unequal weights remember the split.
The paired-switch rank target is

```text
rank(rho_z restricted to ker(mu_-)) >= kappa_car
for generic z.                                        (GRR)
```

Since the target of `rho_z` is the single positive `{a,f}` channel, `(GRR)`
immediately implies `(PSR)`.  The strict C++ lift verifier evaluates eight
deterministic Vandermonde specializations

```text
z_alpha=(mask(alpha)+1)^m,       1<=m<=8.
```

In every one of the 300 cases through label eight, the maximum residue rank
was at least `kappa_car`; hence the strengthened sweep passed.  For example,
at `[3,3,1,1,2,2,2]` the full negative kernel, allocation kernel, carrier
kernel, weighted residue rank, and reservoir dimension are respectively

```text
7, 2, 5, 6, 7.
```

This identifies the exact determinant which an unbounded proof must show is
not identically zero.  The computation does not prove `(GRR)`.

In fact no search over weights was needed in the tested range.  Let `C` be
the unweighted split-forgetting `{a,f}` contraction on `O_-`, and let `D`
be diagonal on the factorization-source summands, multiplying the summand
whose three-element block has bit mask `m` by `m+1`.  Then the first
Vandermonde specialization is simply

```text
rho_1=C D.
```

All 300 cases already satisfy

```text
rank(C D restricted to ker(mu_-)) >= kappa_car.       (LRR)
```

Since `kappa_car=rank(iota_-)-rank(mu_-)`, the standard rank formula for a
stacked map makes `(LRR)` equivalent to the single fixed inequality

```text
rank(mu_- direct_sum C D) >= rank(iota_-).             (P12F.35)
```

Thus the seven-factor doubled-fundamental problem has been reduced to a
rank comparison between two completely explicit maps.  The other seven
powers were retained only as an audit; in every displayed case all eight
had the same rank.  A proof of `(P12F.35)`, for example by a leading-term
degeneration from the marked-edge map `iota_-`, would prove `(PSR)`.

The mask weight can be reduced further.  Split `(P12F.32)` as

```text
O_-=O_q direct_sum O_r
```

according to the odd terminal in the three-element block, and let `C_q`
and `C_r` be the `{a,f}` residual contraction restricted to those two
summands.  The unweighted contraction is `C_q+C_r`; because it factors
through `mu_-`, on `ker(mu_-)` one has

```text
C_q=-C_r.                                             (P12F.36)
```

Weighting only by the minus terminal therefore has, up to a nonzero scalar,
the same restriction as `C_q`.  In all 300 cases the still simpler estimate

```text
rank(C_q restricted to ker(mu_-)) >= kappa_car        (OSR)
```

passed.  Equivalently,

```text
rank(mu_- direct_sum C_q) >= rank(iota_-).             (P12F.37)
```

There is an additional one-sided structural fact, but it has an essential
bottom-label qualification.  Write `mu_s` for the restriction of `mu_-` to
`O_s`, and write `u` for the other odd terminal.  Since `C_s` is ordinary
contraction of the resulting invariant, it factors through `mu_s`;
augmenting by `C_s` cannot increase rank.

**Lemma 12F7 (cut profile and raw leading monomial).**  Order the seven
vertices as

```text
s,a,x,y,z,u,f,
```

where `x<=y<=z` are the three even labels.  For a bracket graph `G`, let
`h_m(G)` be the number of edges crossing the cut after vertex `m`.  In the
coordinate term order whose bracket leading term is

```text
in(p_(ij))=X_i Y_j,                 i<j,
```

the exponent of `X_m` in `in(p_G)` is

```text
(h_m(G)-h_(m-1)(G)+d_m)/2,          h_(-1)=0,           (P12F.39)
```

where `d_m` is the degree at vertex `m`.  Consequently two graphs of the
same multidegree have distinct cut profiles if and only if their raw leading
coordinate monomials are distinct.

**Proof.**  If `o_m` and `i_m` are respectively the numbers of edges from
`m` to later and earlier vertices, then `d_m=o_m+i_m` and
`h_m-h_(m-1)=o_m-i_m`.  Solving gives `(P12F.39)`, while the chosen bracket
term contributes `X_m` exactly on the `o_m` later edges.  QED.

Choose on each four-leg complement its noncrossing bracket basis in the
induced order.  Such a basis vector is indexed by the label on its middle
cut.  Denote by `A_i(lambda)` the source whose triple contains `a`, and by
`F_i(nu)` the source whose triple contains `f`.  Directly adding the unique
three-leg profile to the four-leg profile gives

```text
A_2(lambda): (s,x,0,y,lambda,1),
A_3(lambda): (s,y,x+y,x,lambda,1),
A_4(lambda): (s,z,x+z,z+lambda,lambda,1),

F_2(nu):     (s,s+1,2,nu+1,u+1,1),
F_3(nu):     (s,s+1,s+nu,nu+1,u+1,1),
F_4(nu):     (s,s+1,s+nu,s+u,u+1,1).                  (P12F.40)
```

Here an `A_i` or `F_i` row occurs only when its omitted even label is
`s-1` or `s+1`; also

```text
lambda in {u-1,u+1},
nu in {e_j-1,e_j+1},
```

with `e_j` the first remaining even vertex in the `F_i` complement.

**Proposition 12F8 (interior one-sided injectivity).**  If `s>=3` and
`u>=3`, then

```text
mu_s:O_s -> M is injective.                              (OSI_int)
```

**Proof.**  First compare the profiles in `(P12F.40)`.  Within the `A`
rows, the second coordinate separates unequal omitted labels.  If the
labels agree, the third coordinate separates `A_2` from the other rows,
and the fourth separates `A_3` from `A_4`; here `lambda>=2`.  Within the
`F` rows, the third coordinate separates `F_2`.  A collision between
`F_3(nu_3)` and `F_4(nu_4)` requires

```text
nu_3=nu_4=nu,             nu+1=s+u.                    (P12F.41)
```

Since `nu<=x+1<=s+2` and `s+u-1>=s+2`, `(P12F.41)` forces

```text
u=3,             x=y=z=s+1,             nu=s+2.        (P12F.42)
```

An `A`--`F` collision would first require the omitted `A` label to be
`s+1`.  The third coordinate rules out `A_2`; for `A_3` it forces
`nu=x+1`, after which the fourth coordinates differ; for `A_4`, comparison
of the fifth coordinate and then the fourth gives a difference of two.
Thus `(P12F.42)` is the only repeated profile.

Outside `(P12F.42)`, Lemma 12F7 gives a different raw leading monomial for
every chosen source-basis vector, proving independence.  In `(P12F.42)` all
three even labels are on the upper branch.  Regard the coordinates of the
`s`-vertex as variables over the rational-function field in all other
coordinates, and put `L_i=p_(s i)`.  A relation has the form

```text
sum_(i=1)^3 L_i^s [p_(a i)P_i+p_(f i)Q_i]=0.           (P12F.43)
```

The three binary forms `L_1^s,L_2^s,L_3^s` are independent for `s>=2`:
after writing `L_i=X-t_iY`, the minor in the coefficients of
`X^s,X^(s-1)Y,X^(s-2)Y` is a nonzero scalar times the Vandermonde product
`product_(i<j)(t_i-t_j)`.  Hence every bracketed coefficient in `(P12F.43)`
vanishes.  The irreducibles `p_(a i)` and `p_(f i)` are coprime.  Therefore
`p_(a i)` divides `Q_i`; but `Q_i` has degree zero at vertex `i`, so
`Q_i=0`, and similarly `P_i=0`.  Each four-leg complement basis is
independent, completing the proof.  QED.

The strict C++ verifier now checks the unaugmented one-sided ranks and the
raw-leading profiles separately.  Both one-sided defects were zero in all
300 swept cases with `q,r>=3`, including repeated even labels, and the
unique profile collision `(P12F.42)` was reproduced at
`[3,3,1,1,4,4,4]`; its global one-sided ranks are still `12/12` as the
Vandermonde argument predicts.

The qualification in Proposition 12F8 cannot be deleted.  Direct exact
ranks give

```text
[q,r,a,f,E]=[1,3,1,1,2,2,2]:     q-side=11/12, r-side=11/12;
[q,r,a,f,E]=[1,1,1,1,2,2,2]:     q-side=10/12, r-side=10/12.            (P12F.44)
```

Thus one-sided injectivity is false on the unreduced algebraic family at a
bottom terminal.  Corollary 12F3A shows why these kernels require no
allocation payment in the theorem: every such family contains the forbidden
plus/minus overlap at label one and has already been reduced to fewer
factors.  The interior hypothesis in Proposition 12F8 is therefore automatic
on the support-disjoint two-odd edge target.

Assume now `q,r>=3`, so Proposition 12F8 applies to both sides, and put

```text
B_s=mu_s(ker(C_s)).
```

If `(x,y)` lies in `ker(mu_- direct_sum C_q)`, then `(P12F.36)` also gives
`C_r y=0`, while `mu_q x=-mu_r y`.  One-sided injectivity makes `x` and `y`
unique.  Hence

```text
ker(mu_- direct_sum C_q) is isomorphic to B_q intersect B_r.          (P12F.38)
```

Since `rank(iota_-)=dim O_--kappa_alloc`, `(P12F.37)` is now equivalent to

```text
dim(B_q intersect B_r) <= kappa_alloc.                (IAB)
```

Thus the carrier problem has not disappeared, but it has been converted
into an intersection bounded by the already isolated allocation kernel.
Proposition 12F8 therefore closes the one-sided part of this reduction in
the interior `q,r>=3`.  An injection from the intersection in `(IAB)` into
the allocation relations would establish the interior seven-factor
doubled-fundamental bound.  By Corollary 12F3A there is no additional
`min(q,r)=1` branch in the support-disjoint induction.

The two equal fundamental vertices give a stronger formulation of the
interior target.  Let `tau` exchange `a` and `f`, and let `O_-^+` be its
fixed subspace.  No negative source graph contains the edge `p_(af)`, so
`tau` acts without an orientation sign on the paired source space.

**Lemma 12F9 (the stacked kernel is the symmetric kernel).**  If `q,r>=3`,
then

```text
ker(mu_- direct_sum C_q)=ker(mu_-) intersect O_-^+.     (P12F.45)
```

**Proof.**  Write a stacked-kernel vector as `(x,y)` in `O_q direct_sum
O_r`.  Equation `(P12F.36)` gives `C_r y=0` from `C_q x=0`.  The kernel of
the alternating contraction `V_1 tensor V_1 -> V_0` is the symmetric
summand `V_2`, so both `mu_q x` and `mu_r y` are fixed by `tau`.
Proposition 12F8 applied separately to `O_q` and `O_r` then gives
`tau x=x` and `tau y=y`.  This proves the forward inclusion.  Conversely,
`tau` preserves the direct sum `O_q direct_sum O_r`; hence a fixed vector
has fixed `q`- and `r`-components.  Their images lie in the symmetric
summand and are killed by `C_q` and `C_r`, proving the reverse inclusion.
QED.

Consequently the following kernel inclusion is stronger than `(IAB)` and
would close it immediately:

```text
ker(mu_-) intersect O_-^+  subset ker(iota_-).          (EAL)
```

Indeed `(P12F.45)` and `(EAL)` give
`ker(mu_- direct_sum C_q) subset ker(iota_-)`, and taking dimensions gives
`(IAB)`.  Notice that this is not the false weighted inclusion discussed
below: that counterexample uses `C D`, which does not project to the
`tau`-fixed channel.

The strict C++ verifier now constructs the source involution, the complete
stacked kernel, both `tau` projections, and their images under `iota_-`.
All 300 cases through label eight with `q,r>=3` satisfied `(EAL)`.  Among
the 23 cases with nonzero carrier quotient, the stacked kernel was nonzero
only for

```text
[q,r,a,f,E]=[3,3,1,1,2,2,2],
```

where it was one-dimensional and its `iota_-` image was zero.  The extracted
twelve-term vector is `tau`-fixed and is a sum of equal-`q,r`
transposition/allocation relations over the three even positions.  This
finite classification suggests the precise proof mechanism for `(EAL)`:
separate the `tau`-fixed source jets unless both odd terminals are three
and all three even labels are two, then use that explicit equal-terminal
allocation relation.  The unbounded separation statement is not yet
proved.

There is a six-vertex formulation which removes the source involution from
the remaining proof.  Identify `a=f=v`.  If `G` is the representative of a
`tau`-orbit whose triple contains `a`, merge the `a`- and `f`-vertices of
`G` and call the resulting degree-two-at-`v` graph `H`.  Let `W` be the free
space on these merged graphs and define

```text
m:W -> R_6,                         m([H])=p_H,

delta([H])
 =sum_(j!=v) h_(vj) (p_H/p_(vj)) e_j.                 (P12F.46)
```

Bracket-orientation signs may be absorbed into the basis vectors `e_j`.

**Lemma 12F10 (polarized six-vertex form of EAL).**  On the interior
`q,r>=3`, diagonal restriction identifies the `tau`-fixed negative source
map with `m`, and identifies `iota_-` on that source with `delta`.  Hence

```text
(EAL)  iff  ker(m)=ker(delta).                          (P12F.47)
```

**Proof.**  A symmetric polynomial of degree one at each of `a,f` is
determined by its diagonal value at `a=f=v`; this is ordinary polarization
in characteristic zero.  If the two fundamental edges of `G` end at
`i,j`, its orbit sum restricts to `2p_(vi)p_(vj)` times the common remaining
coefficient.  Stripping the `f`-edge before restriction gives

```text
p_(vi)e_j+p_(vj)e_i,
```

with the evident factor two when `i=j`.  This is exactly `delta([H])`.
Thus the two pairs of maps are identified.  Moreover the column map `phi`
satisfies `phi delta=2m`, so `ker(delta)` is always contained in `ker(m)`;
the reverse inclusion is precisely `(EAL)`.  QED.

The focused C++ mode `--tau-even-sweep` constructs only `(P12F.46)`, without
the much larger `F_0,F_1` carrier complexes.  The strict one-thread runs

```text
maximum label 12: 840 cases;
maximum label 16: 3,360 cases
```

both passed `(P12F.47)`.  In the entire latter domain the only nonzero
kernel of either map was

```text
[q,r,a,f,E]=[3,3,1,1,2,2,2],       dim ker(m)=dim ker(delta)=1.        (P12F.48)
```

Thus the interior problem is now the following exact unbounded statement:
every relation among the structured merged monomials in `(P12F.46)` is
annihilated by the formal `v`-edge derivative.  Equivalently, those
relations use only Plucker moves not involving the merged vertex, except
for the explicit allocation relation `(P12F.48)`.

One large subfamily admits a genuine divisor-separation certificate.  For a
finite family of bracket graphs, choose a bracket `p_(ij)`, group the graphs
by its exponent, divide each group by that common power, and then impose
`p_(ij)=0`, which merges vertices `i,j`.  Repeat on each group.

**Lemma 12F11 (recursive divisor separation).**  If this recursion reaches
families of distinct graphs on at most three vertices, then the original
bracket monomials are linearly independent.

**Proof.**  In a putative relation, let `m` be the smallest exponent of
`p_(ij)`.  Divide by `p_(ij)^m` and restrict to the irreducible divisor
`p_(ij)=0`.  All terms with larger exponent vanish; the terms with exponent
`m` give a relation among the merged graphs.  Recursive independence kills
their coefficients.  Remove them and repeat with the next exponent.  On at
most three vertices the bracket ring is a polynomial ring, so distinct
bracket monomials are independent.  QED.

The focused verifier constructs this recursion exactly.  The full rank sweep
through label sixteen certified 3,325 of the 3,360 cases directly.  Of the 35
unresolved by this sufficient criterion, 34 had `q=r`.  A faster
unequal-terminal-only divisor sweep through label forty checked 263,340
parameter tuples, of which 75,972 had nonempty source, and again found the
sole unequal-terminal exception

```text
[q,r,a,f,E]=[3,5,1,1,2,4,4],
```

and a separate circular standard-monomial order gives ten distinct
leading terms there.  In the merged vertex order `(q,r,v,e_1,e_2,e_3)`, the
successful circular order is

```text
(q,e_1,e_2,r,e_3,v).                                  (P12F.49)
```

Thus every unequal-terminal case in that domain has an explicit algebraic
independence certificate, rather than only a modular rank computation.  The
following formulas reduce the unbounded statement to a small
standard-profile table.

Write the merged vertices as `(q,r,v,e_1,e_2,e_3)`.  Let `s` be one of
`q,r`, let `u` be the other terminal, choose `z=e_i=s-1` or `s+1`, and write
the two remaining even labels in their induced order as `x<=y`.  Put

```text
A_s^-(z)=p_(sv)p_(sz)^(s-1),              z=s-1,
A_s^+(z)=p_(sz)^s p_(vz),                 z=s+1.       (P12F.50)
```

For the two possible complement channels put

```text
B_u^-(x,y)
 =p_(uv)
  p_(ux)^((u-1+x-y)/2)
  p_(uy)^((u-1-x+y)/2)
  p_(xy)^((x+y-u+1)/2),                              (P12F.51)
```

and, when `y-x<=u-1`,

```text
B_u^+(x,y)
 =p_(ux)^((u+x-y-1)/2)
  p_(uy)^((u-x+y+1)/2)
  p_(vx)
  p_(xy)^((x+y-u-1)/2).                              (P12F.52)
```

At the remaining upper wall `y-x=u+1`, replace `(P12F.52)` by

```text
B_u^+(x,y)=p_(uy)^u p_(vy)p_(xy)^x.                  (P12F.53)
```

A displayed monomial is included only when all its exponents are
nonnegative.

**Lemma 12F12 (explicit merged source list).**  The merged graph family in
`(P12F.46)` consists exactly, with every factorization-source index retained,
of

```text
A_s^epsilon(z) B_u^eta(x,y),
s in {q,r},  z=e_i=s+epsilon,  eta in {-,+},           (P12F.54)
```

whenever the indicated complement channel is admissible.

**Proof.**  The three-leg invariant on `(s,v,z)` is multiplicity-free.
Solving its three degree equations gives `(P12F.50)`.  On the ordered
four-leg complement `(u,v,x,y)`, let `t` be the label crossing the middle
cut.  Since `v` has degree one, `t=u-1` or `u+1`.  For `t=u-1`, solving the
four vertex-degree equations gives `(P12F.51)`.  For `t=u+1`, noncrossing
forces the unique `v`-edge to end at `x` unless `y-x=u+1`; these two cases
give `(P12F.52)` and `(P12F.53)`.  Conversely each nonnegative displayed
exponent vector has the required degrees and is noncrossing on the
four-leg complement, so it is one of its standard basis vectors.  QED.

For a circular vertex order `O`, straighten a bracket monomial into the
noncrossing standard basis and denote the least and greatest monomials in
its support by `ell_O` and `g_O`.

**Lemma 12F13 (extremal standard-monomial criterion).**  A finite family of
bracket monomials is linearly independent if either its `ell_O` values or
its `g_O` values are pairwise distinct for one order `O`.

**Proof.**  In a putative relation choose the largest source with respect to
the relevant distinct extremal standard monomial.  That basis monomial
occurs with a positive integer coefficient in its own Plucker
straightening and in no source with a smaller extremum, so its coefficient
in the relation cannot cancel.  Descending induction kills every source
coefficient.  The argument for the least extremum is identical.  QED.

**Lemma 12F13A (nested cut extremum).**  Fix a vertex order and let `h_j(G)`
be the number of edges of `G` crossing the prefix cut after position `j`.
The positive Plucker straightening of `p_G` contains a unique noncrossing
graph `N(G)` with

```text
h_j(N(G))=h_j(G)                 for every j.          (P12F.54a)
```

Every other graph `H` in its support satisfies `h_j(H)<=h_j(G)` for every
`j`, with a strict inequality for at least one `j`.  Consequently a family
with pairwise distinct cut profiles is linearly independent.

**Proof.**  For `a<b<c<d`, the nested resolution

```text
p_(ac)p_(bd) -> p_(ad)p_(bc)
```

preserves every prefix cut weight, whereas the separated resolution

```text
p_(ac)p_(bd) -> p_(ab)p_(cd)
```

lowers by two precisely the cuts from `b` through `c-1`.  Repeated nested
resolution terminates: the sum of the squared edge lengths strictly
increases and is bounded in the fixed multidegree.  A noncrossing graph is
uniquely determined by its degrees and cut profile.  Indeed `(P12F.39)`
determines the number of right-going half-edges at every vertex, and
noncrossing forces the open half-edges to close in last-in-first-out order.
This proves existence and uniqueness of `N(G)` and the componentwise bound
for every other resolution branch.

In a relation among sources with distinct cut profiles, choose a source
whose profile is maximal in the componentwise partial order.  Its `N(G)`
cannot occur in the straightening of any other source and has positive
coefficient, so that source coefficient vanishes.  Iteration proves
independence.  QED.

**Lemma 12F14 (lower-wall triple family).**  If

```text
r=2q-1,                 e_1=e_2=e_3=q-1,              (P12F.54b)
```

then the unequal-terminal merged source monomials are linearly independent.

**Proof.**  Only a `q`-side source is possible.  Its selected three-leg
channel is the lower neighbor `q-1`.  On the complement, the two even labels
are both `q-1`; their maximum fusion channel is `2q-2=r-1`, so only the lower
complement channel occurs.  Thus there are exactly three sources,

```text
G_i=p_(qv)p_(rv) p_(q e_i)^(q-1)
       product_(j!=i) p_(r e_j)^(q-1),                (P12F.54c)
```

one for each even vertex.  Regard the coordinates of the `q`-vertex as
variables over the rational-function field in all other coordinates.  After
removing the common factor `p_(qv)p_(rv)`, a relation among the `G_i` is a
relation among the three powers `p_(q e_i)^(q-1)` with nonzero coefficients
in that field.  Since `q>=3`, these are three powers of distinct generic
binary linear forms of degree at least two.  Their first three coefficient
rows have determinant equal, up to a nonzero scalar, to the Vandermonde
product of the three slopes.  They are therefore independent, and all three
source coefficients vanish.  QED.

**Lemma 12F15 (the unequal cubic chart).**  The merged source monomials for

```text
(q,r;e_1,e_2,e_3)=(3,5;4,4,4)                        (P12F.54d)
```

are linearly independent.

**Proof.**  Use the ordered vertices

```text
(0,1,2,3,4,5)=(q,e_3,r,e_1,e_2,v).
```

In the edge order

```text
01,02,03,04,05,12,13,14,15,23,24,25,34,35,45,
```

the greatest noncrossing monomials in the positive Plucker straightenings
are as follows.  The source row `(s,i,eta)` means that terminal `s` selects
`e_i` and the complement uses channel `eta`.

```text
(q,1,+): (1,2,0,0,0,3,0,0,0,0,0,0,3,1,1)
(q,1,-): (2,1,0,0,0,2,0,0,0,2,0,0,2,0,2)
(q,2,+): (0,0,1,2,0,3,1,0,0,2,0,0,0,0,2)
(q,2,-): (1,0,0,2,0,2,1,0,0,3,0,0,0,0,2)
(q,3,+): (3,0,0,0,0,0,0,0,1,3,2,0,1,0,1)
(q,3,-): (3,0,0,0,0,0,0,0,1,2,2,1,2,0,0)
(r,1,+): (3,0,0,0,0,0,0,1,0,4,1,0,0,0,2)
(r,1,-): (2,0,0,0,1,0,0,2,0,4,1,0,0,0,1)
(r,2,+): (3,0,0,0,0,1,0,0,0,3,1,0,1,0,2)
(r,2,-): (2,0,0,0,1,2,0,0,0,2,1,0,2,0,1)
(r,3,+): (0,1,2,0,0,4,0,0,0,0,0,0,2,0,2)
(r,3,-): (0,1,1,0,1,4,0,0,0,0,0,0,3,0,1).       (P12F.54e)
```

To verify the table, substitute `(P12F.50)`--`(P12F.53)` and recursively use

```text
p_(ac)p_(bd)=p_(ab)p_(cd)+p_(ad)p_(bc),   a<b<c<d.
```

At each recursion node retain the lexicographically greatest of the two
already-straightened children.  The standard Plucker straightening used in
Proposition 12F4 is a terminating recurrence into the noncrossing basis;
evaluating that recurrence on these twelve fixed inputs gives exactly the
vectors in `(P12F.54e)`.  They are pairwise distinct, so Lemma 12F13 proves
independence.  QED.

The remaining unequal-terminal step is now the following explicit
assumption; it has not yet been proved for unbounded labels.

**Assumption 12F16 (two-chart profile separation).**  Let `3<=q<r` be odd
and let `e_1<=e_2<=e_3` be even.  For the admissible monomials `(P12F.54)`,
excluding the families in Lemmas 12F14 and 12F15, `g_O` is injective in the
chart prescribed below:

```text
O_1=(q,e_1,e_2,r,e_3,v),
  unless one of the next two lines applies;

O_2=(q,r,e_3,e_1,v,e_2),
  if e_1=e_2=q-1;
                                                        (P12F.55)
```

By Lemmas 12F13--12F15, Assumption 12F16 implies that `m` is injective
whenever `q!=r`, and hence proves `(P12F.47)` in the unequal-terminal branch.

The strict C++ chart audit propagates Plucker *supports* as Boolean sets, so
this check has neither modular arithmetic nor coefficient-overflow risk.  It
tested the two charts in `(P12F.55)`, together with the cubic chart of Lemma
12F15, through label twenty-four: 9,190 nonempty unequal-terminal families,
no uncovered family and no partition failure.  The first chart handled
8,914, the second 275, and the cubic chart only `(3,5;4,4,4)`.  The only
families through label sixteen in which the least, rather than greatest,
extremum was required were exactly the lower-wall triples `(P12F.54b)`, now
covered by Lemma 12F14.  This is finite evidence for Assumption 12F16; the
required proof is the unbounded comparison of the explicit exponent profiles
`(P12F.50)`--`(P12F.53)` under the two orders `(P12F.55)`.  A simpler proposed
shortcut is false: always choosing the separated child at every crossing did
not produce `g_O` in 44 of the 435 nonempty families through label twelve, so
the proof must retain the genuine extremal straightening recurrence.

There is a second, substantially smaller, route to the same unequal-terminal
injectivity.  Work over the fraction field `K` of the bracket algebra on the
four nonterminal vertices `v,e_1,e_2,e_3`.  No source in `(P12F.54)` contains
`p_(qr)`, so it has a unique factorization

```text
p_G=C_G F_G(q) H_G(r),                              (P12F.56)
```

where `C_G` is a nonzero element of `K`, while `F_G` and `H_G` are binary
forms of degrees `q` and `r`, respectively.

**Lemma 12F17 (two-terminal specialization criterion).**  If the pure
tensors

```text
T_G=F_G tensor H_G in Sym^q(K^2) tensor Sym^r(K^2)   (P12F.57)
```

are linearly independent over `K`, then the merged source monomials `p_G`
are linearly independent over the ground field.

**Proof.**  A ground-field relation among the `p_G` remains a relation after
extending scalars to `K`.  By `(P12F.56)` its rows are the nonzero scalar
multiples `C_G T_G`.  Independence of the `T_G` therefore kills every
coefficient.  QED.

The terminal-tensor test has one transparent exceptional family, which can
be closed without straightening.

**Lemma 12F18 (the terminal-tensor circuit family).**  For every odd
`3<=q<r`, the merged source monomials with

```text
(e_1,e_2,e_3)=(2,q+1,r+1)                           (P12F.58)
```

are linearly independent.

**Proof.**  Put `c=e_1`, `a=e_2`, and `b=e_3`.  Substitution in
`(P12F.50)`--`(P12F.53)` always gives the following four sources:

```text
M_1=p_(qa)^q p_(rb)^r p_(va)p_(vc)p_(cb),
M_2=p_(qa)^q p_(rv)p_(rb)^(r-1)p_(va)p_(cb)^2,
M_3=p_(qa)^q p_(rb)^r p_(vb)p_(vc)p_(ca),
M_4=p_(qv)p_(qa)^(q-1)p_(rb)^r p_(vb)p_(ca)^2.       (P12F.59)
```

The four-source subfamily is independent.  Indeed, after removing the common terminal factor
`p_(qa)^(q-1)p_(rb)^(r-1)`, their terminal tensors are, in order,

```text
p_(qa)p_(rb),  p_(qa)p_(rv),
p_(qa)p_(rb),  p_(qv)p_(rb).                         (P12F.60)
```

The three distinct tensors in `(P12F.60)` are independent because
`p_(qa),p_(qv)` are independent linear forms in `q` and
`p_(rb),p_(rv)` are independent linear forms in `r`.  Thus a relation among
these four first kills the coefficients of `M_2` and `M_4`.  The remaining
coefficient is

```text
c_1 p_(va)p_(vc)p_(cb)+c_3 p_(vb)p_(vc)p_(ca).       (P12F.61)
```

The two bracket monomials are not constant multiples: for example the
irreducible factor `p_(va)` divides only the first.  Hence `c_1=c_3=0`.

If `r=q+2`, the lower `r`-channel at `a=r-1` gives one additional source,

```text
M_5=p_(qb)^q p_(rv)p_(ra)^(r-1)p_(vb)p_(cb)^2.       (P12F.62)
```

If `q=3`, the lower `q`-channel at `c=q-1` also gives

```text
M_6=p_(qv)p_(qc)^2 p_(ra)p_(rb)^(r-1)p_(va)p_(ab)^2,
M_7=p_(qv)p_(qc)^2 p_(rv)p_(ra)p_(rb)^(r-2)p_(ab)^3. (P12F.63)
```

There are no other sources.  To remove these additions, first set `q`
proportional to `a`; the four sources `(P12F.59)` vanish.  If `q>=5`, only
`M_5` can remain, so its coefficient vanishes.  If `q=3` and `r>=7`, only
`M_6,M_7` remain; after removing their common `r`-factor
`p_(ra)p_(rb)^(r-2)`, their residual linear forms are `p_(rb)` and
`p_(rv)`, so both coefficients vanish.

Finally suppose `q=3,r=5`.  All three additions remain, and their `r`-factors
are

```text
p_(rv)p_(ra)^4,  p_(ra)p_(rb)^4,
p_(rv)p_(ra)p_(rb)^3.                                (P12F.64)
```

Setting `r` proportional to `b` first isolates `M_5`; after it is removed,
the other two factors in `(P12F.64)` have the common factor
`p_(ra)p_(rb)^3` and the independent residual linear forms `p_(rb)` and
`p_(rv)`.  Thus every additional coefficient vanishes in all cases.  The
four-source argument `(P12F.60)`--`(P12F.61)` finishes the proof.  QED.

This leaves the following alternative to Assumption 12F16.

**Assumption 12F19 (two-terminal tensor separation).**  Outside the family
`(P12F.58)`, the tensors `(P12F.57)` attached to the admissible sources
`(P12F.54)` are linearly independent over `K`.

Lemma 12F17 and Lemma 12F18 show that Assumption 12F19 also proves
injectivity of `m` for every `q!=r`.  It is a smaller interpolation problem
than Assumption 12F16: only two binary variables remain, and every row is a
pure tensor of products of four fixed linear forms.

One infinite family missed by the simplest one-root certificate is also
elementary.

**Lemma 12F20 (adjacent double-upper family).**  Let `q>=5`, let `r=q+2`,
let `a_1,a_2` both have label `q+1`, and let `c` have label `2q`.  Then the
merged sources are linearly independent.

**Proof.**  Ignore nonzero scalars in `K`.  The four sources in which `q`
selects one of the two vertices `a_i` have terminal tensors

```text
p_(q a_i)^q tensor
 p_(r a_j)p_(rc)^q {p_(rc),p_(rv)},       {i,j}={1,2}. (P12F.65)
```

The four sources in which the lower `r`-channel selects `a_i` have terminal
tensors

```text
p_(qc)^(q-1){p_(qc),p_(qv)} tensor
 p_(rv)p_(r a_i)^(q+1),                    i in {1,2}. (P12F.66)
```

There are no other sources.  Set `r` proportional to `c`.  Every tensor in
`(P12F.65)` vanishes, while those in `(P12F.66)` survive.  The two `r`-forms
`p_(rv)p_(r a_i)^(q+1)` are independent, and the two residual `q`-forms
`p_(qc),p_(qv)` are independent, so the four coefficients in `(P12F.66)`
vanish.  In `(P12F.65)`, the powers `p_(q a_1)^q,p_(q a_2)^q` are
independent; for each fixed `i`, the two residual `r`-forms
`p_(rc),p_(rv)` are independent.  Thus its four coefficients vanish as
well.  Lemma 12F17 finishes the proof.  QED.

The jet certificate used above has a short abstract formulation.

**Lemma 12F21 (one-root jet separation).**  Fix a nonterminal vertex `z`.
For each integer `h`, group the tensors `(P12F.57)` for which `F_G(q)` has
vanishing order `h` on `p_(qz)=0`.  If, in every group, the forms `H_G(r)`
are independent, then all the tensors are independent.  The same statement
holds with `q` and `r` exchanged.

**Proof.**  Use `p_(qz)` as a local parameter at a generic point of the
divisor.  In a relation, take the smallest vanishing order `h`.  The
coefficient of that order is a relation among the corresponding `H_G`,
multiplied by nonzero elements of `K`; hence every coefficient in the group
vanishes.  Iterate through the vanishing orders.  QED.

The remaining fixed jet exceptions also have short terminal proofs.

**Lemma 12F22 (adjacent repeated-neighbor families).**  Let `r=q+2`.  The
merged sources are linearly independent in either of the following families:

```text
(e_1,e_2,e_3)=(q+1,q+1,q+1), q>=5;
(e_1,e_2,e_3)=(q+1,q+1,q+3), q>=3.                  (P12F.67)
```

**Proof.**  First take the all-repeated family, call the three even vertices
`a,b,c`, and put `m=(q+1)/2>=3`.  Setting `r` proportional to `v` isolates
the three upper
`q`-channel sources that do not contain `p_(rv)`.  Their `q`-factors are
`p_(qa)^q,p_(qb)^q,p_(qc)^q`, which are independent by the Vandermonde
argument of Lemma 12F14, so their coefficients vanish.  Divide the remaining
relation by its common factor `p_(rv)`.  Its six distinct `r`-factors are

```text
p_(ra)^(2m),p_(rb)^(2m),p_(rc)^(2m),
p_(ra)^m p_(rb)^m,p_(ra)^m p_(rc)^m,p_(rb)^m p_(rc)^m. (P12F.68)
```

These six binary forms are independent.  Indeed, send the three linear forms
projectively to `X,Y,X+Y`.  For a relation among

```text
X^(2m),Y^(2m),(X+Y)^(2m),X^mY^m,
X^m(X+Y)^m,Y^m(X+Y)^m,
```

the coefficients of `X^(2m-1)Y` and `X^(2m-2)Y^2` kill the
`(X+Y)^(2m)` and `X^m(X+Y)^m` terms: the resulting two-by-two determinant is
nonzero (eliminating the first row leaves the factor `m^2`).  The symmetric
two coefficients near the other end kill `Y^m(X+Y)^m`; the middle coefficient and the two endpoint
coefficients kill the remaining three terms.  For each pure power in
`(P12F.68)`, its two accompanying `q`-forms are independent; each mixed
factor accompanies only one remaining source.  Hence all coefficients
vanish.

Now take the second family and call the repeated vertices `a,b` and the last
vertex `c`.  Setting `r` proportional to `c` kills the four `q`-selected
sources and the two upper-`r` sources.  The four surviving lower-`r` sources
have independent `r`-factors
`p_(rv)p_(ra)^(q+1),p_(rv)p_(rb)^(q+1)`; over each factor their two residual
`q`-forms are `p_(qc),p_(qv)`, so all four coefficients vanish.  In the
remaining relation, set `q` proportional to `a`.  This isolates the two
sources in which `q` selects `b`; their residual `r`-forms are
`p_(rc),p_(rv)`, so both coefficients vanish.  Setting `q` proportional to
`b` similarly kills the two sources selecting `a`.  The last two sources
share their powers of `p_(rc),p_(qa),p_(qb)` and have independent residual `q`-forms
`p_(qb),p_(qv)`.  Their coefficients vanish.  Lemma 12F17 applies in both
families.  QED.

Consequently the alternative unequal-terminal bottleneck is only the
following finite-pattern statement.

**Assumption 12F23 (jet-classification statement).**  Outside

```text
(2,q+1,r+1),
(2,q+1,q+1) with r=q+2,
(q+1,q+1,t) with r=q+2 and t in {q+1,q+3,2q},
(3,5;4,4,4),                                         (P12F.69)
```

at least one of the eight applications of Lemma 12F21 obtained by choosing
the terminal `q` or `r` and one of the four roots `v,e_1,e_2,e_3` separates
the tensors `(P12F.57)`.  This statement depends only on comparing the
explicit root multiplicities in `(P12F.50)`--`(P12F.53)`.  Lemmas 12F15,
12F17--12F22, and 12F24 show that Assumption 12F23 proves the entire
unequal-terminal branch.

**Lemma 12F24 (adjacent bottom-double-upper family).**  Let `r=q+2`, let
`c` have label two, and let `a_1,a_2` both have label `q+1`.  Then the merged
sources are linearly independent.

**Proof.**  Suppose first that `q>=5`.  Up to nonzero factors in `K`, the
four sources in which `q` selects `a_i` have terminal tensors

```text
p_(q a_i)^q tensor
 p_(rc)p_(r a_j)^q {p_(rv),p_(r a_j)},    {i,j}={1,2}. (P12F.68a)
```

The four sources in which the lower `r`-channel selects `a_i` have terminal
tensors

```text
p_(q a_j)^(q-1){p_(qv),p_(q a_j)} tensor
 p_(rv)p_(r a_i)^(q+1),                   {i,j}={1,2}. (P12F.68b)
```

Set `r` proportional to `c`.  The first four tensors vanish and the last
four survive.  Their two `r`-factors are independent powers, and over each
one the residual `q`-forms `p_(qv),p_(q a_j)` are independent.  Thus their
coefficients vanish.  The powers `p_(q a_1)^q,p_(q a_2)^q` then separate the
first four tensors, and over each power the residual `r`-forms in
`(P12F.68a)` are independent.

If `q=3`, the lower `q`-channel at `c` contributes two further sources.  They
have common `q`-factor `p_(qv)p_(qc)^2` and, up to nonzero factors in `K`,
their `r`-factors are

```text
p_(rv)p_(r a_1)^2p_(r a_2)^2,
p_(r a_1)^2p_(r a_2)^3.                              (P12F.68c)
```

After setting `r` proportional to `c`, set `q` proportional to `c`.  The
two additional sources vanish.  The four sources `(P12F.68b)` have `q`-forms

```text
p_(qv)p_(q a_1)^2, p_(q a_1)^3,
p_(qv)p_(q a_2)^2, p_(q a_2)^3.                     (P12F.68d)
```

They are independent: send `p_(q a_1),p_(q a_2),p_(qv)` to `X,Y,X+Y` and
compare the coefficients of `X^3,X^2Y,XY^2,Y^3`.  Hence their coefficients
vanish.  The two forms `(P12F.68c)` have common factor
`p_(r a_1)^2p_(r a_2)^2` and independent residual forms
`p_(rv),p_(r a_2)`, so their coefficients also vanish.  The first four
sources are handled as above.  Lemma 12F17 finishes the proof.  QED.

**Lemma 12F24A (the full adjacent two-neighbor family).**  Let `r=q+2`
and suppose that two of the three even vertices have label `q+1`.  The
third even label may be any positive even integer `t`.  Then the merged
sources are linearly independent.

**Proof.**  Call the repeated vertices `a,b` and the third vertex `c`.
The cases `t=2`, `t=q+1`, and `t=q+3` are Lemma 12F24, Lemma 12F22
(together with Lemma 12F15 when `q=3,t=4`), and Lemma 12F22,
respectively.

Suppose next that `t>=4` and that `t` is not one of `q-1,q+1,q+3`.
Then only `a,b` can be selected by either odd terminal.  Every source in
which `q` selects `a` or `b` has its `r`-form divisible by
`p_(rc)^(t/2)`.  The sources in which `r` selects `a` or `b` have
`r`-forms

```text
p_(rv)p_(ra)^(q+1), p_(rv)p_(rb)^(q+1).              (P12F.68e)
```

Their images modulo `p_(rc)^2` are independent: their two-by-two jet
determinant at `p_(rc)=0` is nonzero because their ratio
`(p_(ra)/p_(rb))^(q+1)` is nonconstant.  Thus the two blocks in
`(P12F.68e)` separate.  Within either block the at most two `q`-forms
`B_q^-` and `B_q^+` are independent, since the former contains `p_(qv)`
and the latter does not.  All coefficients of the `r`-selected sources
therefore vanish.  The remaining sources have the independent `q`-powers
`p_(qa)^q,p_(qb)^q`; within either power, `B_r^-` and `B_r^+` are
independent for the same `p_(rv)` divisibility reason.  Hence these
coefficients vanish as well.  The case `t=2` used above covers the only
case in which the two-jet quotient is unavailable.

It remains to take `t=q-1` with `q>=5`.  Put `m=(q+1)/2`, so
`q=2m-1` and `r=2m+1`.  Besides the eight sources just considered, the
lower `q`-channel can select `c`.  Its two sources have the common
`q`-form

```text
F_0=p_(qv)p_(qc)^(2m-2)                              (P12F.68f)
```

and, up to nonzero elements of `K`, the two `r`-forms

```text
H_0=p_(ra)^m p_(rb)^(m+1),
H_1=p_(rv)p_(ra)^m p_(rb)^m.                         (P12F.68g)
```

The four sources in which `r` selects `a` or `b` have, after removing
the common factor `p_(qc)^(m-2)`, the following four `q`-forms:

```text
p_(qv)p_(qb)^m, p_(qb)^(m+1),
p_(qv)p_(qa)^m, p_(qa)^(m+1).                        (P12F.68h)
```

Together with the residual form `p_(qv)p_(qc)^m` from `(P12F.68f)`,
these five forms are independent.  Indeed, send `p_(qa),p_(qb),p_(qv)`
to `X,Y,X+Y` and write `p_(qc)=X+lambda Y`, with `lambda` generic.  The
four forms `(P12F.68h)` span the four monomials with `X`-exponents
`0,1,m,m+1`, whereas

```text
(X+Y)(X+lambda Y)^m
```

has a nonzero generic coefficient at `X^(m-1)Y^2`; here `m>=3`.

Evaluate a putative tensor relation at `r=c`.  The four independent
forms `(P12F.68h)` kill the four `r`-selected coefficients, while the
coefficient of `F_0` gives one relation between the two coefficients in
`(P12F.68g)`.  Differentiate once in the local parameter `p_(rc)` and
evaluate again.  The four sources in which `q` selects `a` or `b` still
vanish because their `r`-forms are divisible by `p_(rc)^(m-1)` and
`m-1>=2`.  The resulting second relation is independent of the first:
the first-jet determinant of `H_0,H_1` is nonzero because
`H_0/H_1=p_(rb)/p_(rv)` is nonconstant.  Hence the two coefficients of
`F_0` also vanish.  The last four sources separate by
`p_(qa)^q,p_(qb)^q` and then by their two nonproportional `r`-forms, as
in the preceding paragraph.  Lemma 12F17 finishes every case.  QED.

**Lemma 12F24B (the cubic-terminal repeated-neighbor families).**  Let
`q=3`.  The merged sources are linearly independent when all three even
labels equal `r-1` or all three equal `r+1`.  They are also independent for
`(q,r;e_1,e_2,e_3)=(3,5;4,4,8)`.

**Proof.**  The case `r=5` with three labels `r-1=4` is Lemma 12F15.
In every other repeated case, only the `r`-terminal can select an even
vertex.  If the selected vertex is `e_i`, the two complement channels have
`q`-forms, up to elements of `K`,

```text
p_(qv)p_(q e_j)p_(q e_k),
p_(q e_j)p_(q e_k)^2,             {i,j,k}={1,2,3},   (P12F.68k)
```

which are nonproportional.  The common `r`-form of this two-source block is

```text
p_(rv)p_(r e_i)^(r-1)    if e_i=r-1,
p_(r e_i)^r              if e_i=r+1.                 (P12F.68l)
```

The three forms on either line of `(P12F.68l)` are independent: after
removing `p_(rv)` in the first line, this is the Vandermonde independence
of three distinct generic powers of degree at least four.  Thus the three
blocks separate, and `(P12F.68k)` separates the two sources in each block.

Finally write the exceptional labels `(3,5;4,4,8)` as `a,b,c`.  The four
sources in which `q` selects `a` or `b` have terminal factors

```text
{p_(qa)^3,p_(qb)^3} tensor {p_(rc)^5,p_(rv)p_(rc)^4}, (P12F.68m)
```

while the two sources in which `r` selects `a` or `b` have terminal factors

```text
p_(qc)^3 tensor {p_(rv)p_(ra)^4,p_(rv)p_(rb)^4}.     (P12F.68n)
```

The four displayed `r`-forms are independent.  Reduce a relation modulo
`p_(rv)` to kill the coefficient of `p_(rc)^5`; after division by
`p_(rv)`, the remaining three fourth powers are Vandermonde-independent.
The `q`-forms within each resulting block are also independent.  Lemma
12F17 completes the proof.  QED.

**Lemma 12F25 (root-merger certificate).**  Consider a family of products
of four generic binary linear forms.  If coalescing some two of the four
roots produces a linearly independent family of products in the remaining
three generic forms, then the original family is linearly independent.

**Proof.**  Any maximal minor of the coefficient matrix is a polynomial in
the four root coordinates.  Independence after coalescence says that one
such polynomial has a nonzero value on that specialization; it therefore is
not the zero polynomial.  QED.

**Lemma 12F25A (the proved `V_2` product criterion).**  Let
`L_1,...,L_n` be generic binary linear forms and let

```text
P_i=product_j L_j^(w_i(j)),        |w_i|=d.           (P12F.68i)
```

Remove the coordinatewise common minimum of all the `w_i`, and relabel the
residual rows and their common degree as `w_i` and `d`.  Suppose that, after
permuting coordinates, every residual vector lies in
`{0,1}^(n-2) times N^2`, and that every nonempty row subset `I` satisfies

```text
|I| + |coordinatewise_min_(i in I) w_i| <= d+1.      (P12F.68j)
```

Then the products `(P12F.68i)` are linearly independent.  The conclusion
also holds for four-root products if some pair of roots can be coalesced
and the merged exponent rows satisfy these conditions.

**Proof.**  After removing the common factor, take `k=d+1` in Theorem 5.1
of Greaves--Syatriadi, *Reed--Solomon codes over small fields with
constrained generator matrices*,
[arXiv:1808.06306](https://arxiv.org/abs/1808.06306).  Each set
`P(k,w_i)` in their notation is a singleton because `|w_i|=k-1`;
their property `V_2(k)` is exactly `(P12F.68j)` together with the displayed
residual-multiplicity condition.  Their theorem therefore gives
independence for generic roots.  The characteristic-zero statement follows
from the same nonzero coefficient determinant (or by reduction at a prime
where that determinant remains nonzero).  The final assertion follows from
Lemma 12F25.  QED.

This gives a still sharper sufficient version of Assumption 12F23.

**Assumption 12F26 (merge-jet classification).**  Outside the families
proved in Lemmas 12F15, 12F18, 12F24A, and 12F24B, at least one of the
eight one-root filtrations of Lemma 12F21 has the following property: in
every vanishing-order group, some pair of the four opposite-terminal roots
can be coalesced so that the resulting three-root products are independent.

By Lemmas 12F17, 12F21, and 12F25, Assumption 12F26 implies Assumption 12F23
and closes the unequal-terminal branch.

The algebraic rank assertion inside Assumption 12F26 can now be replaced by
a published theorem.  The remaining formulation is purely combinatorial.

**Assumption 12F26A (`V_2` jet classification).**  Outside the families in
Lemmas 12F15, 12F18, 12F24A, and 12F24B, one of the eight one-root
filtrations of Lemma 12F21 has the following property in every
vanishing-order group: either the four-root exponent rows themselves, or
the rows obtained by coalescing some pair of roots, satisfy the subset and
residual-multiplicity hypotheses of Lemma 12F25A.

Lemma 12F25A shows directly that Assumption 12F26A implies Assumption 12F26.
Unlike the earlier determinant formulation, 12F26A requires only integer
comparisons among the exponents in `(P12F.50)`--`(P12F.53)`.

The `q=3,r>=7` portion of this classification has an exact unbounded
Presburger certificate.  In half-label variables, partition each ordered
even label into

```text
1, 2, [3,r_half-1], r_half, r_half+1, [r_half+2,infinity].
```

The 56 ordered category triples exhaust the domain.  For each triple the
C++ verifier constructs precisely the active rows of `(P12F.50)`--
`(P12F.53)`, negates all four sufficient one-root `V_2` certificates, and
submits the resulting quantifier-free linear-integer formula to Z3's
`qflia` tactic.  All 56 formulas are unsatisfiable.  The durable transcript
is `certificates/su2_v2_q3_large_chambers.log`; its SHA-256 digest is
`27c3c96ef136e83f5909c2e5f76723972300521d02bb78dba690040652ec7910`.
This is a computer-assisted proof for the entire `q=3,r>=7` branch, not a
finite label sweep.  The affine source encoding was independently compared,
with multiplicities, against the noncrossing graph generator in all 840
unequal- and equal-terminal tuples through label twelve (585 nonempty); there
were no mismatches.  The remaining fixed branch `q=3,r=5` has a separate
20-chamber exact replay: nineteen chambers have a uniform one-root `V_2`
certificate and the last is empty.  Its durable transcript is
`certificates/su2_v2_q3_r5_chambers.log`, with hashes in the companion
metadata.  Thus the entire `q=3` unequal-terminal branch is discharged.

The adjacent branch `r=q+2`, `q>=5`, also has an exact unbounded chamber
certificate.  Partition each ordered even label into

```text
<=q-1, q+1, r, r+1, >=q+3
```

in half-label coordinates.  The resulting 35 ordered triples exhaust the
domain.  Of these, 26 have a single uniform root-filtration certificate,
four are empty, and five are precisely the adjacent two-neighbor family of
Lemma 12F24A.  Hence all 35 chambers are closed.  The durable transcript is
`certificates/su2_v2_adjacent_chambers.log`; its companion metadata records
the executable hash.  This discharges the entire adjacent branch of
Assumption 12F26A.

Finally, the separated branch `q>=5`, `r>=q+4` has an exact unbounded
certificate.  Put `Q=(q-1)/2`, `R=(r-1)/2`; each ordered even half-label lies
in

```text
<=Q-1, Q, Q+1, [Q+2,R-1], R, R+1, >=R+2.
```

The seven categories give 84 ordered triples.  The verifier finds 74 uniform
one-root `V_2` certificates and ten empty chambers.  The durable transcript
is `certificates/su2_v2_separated_chambers.log`, with source, executable, and
log hashes in its metadata.  Together with the `q=3` and adjacent results,
this proves Assumption 12F26A and therefore closes the entire unequal-terminal
part of `(P12F.47)`.

The strict C++ terminal-factor audit specializes the four fixed linear forms
at distinct integers and computes the resulting pure-tensor rank over the
prime field of order `1,000,000,007`.  A full rank is an exact certificate:
a nonzero minor modulo the prime is a nonzero integral minor.  Through label
forty, 75,801 of the 75,972 nonempty unequal-terminal families had full
terminal-tensor rank.  The remaining 171 were *exactly* `(P12F.58)`, with a
one-dimensional terminal-tensor circuit, and are now proved independently by
Lemma 12F18.  The one-root test of Lemma 12F21 certified 75,789 families.
Its other twelve failures were exactly eight instances of Lemma 12F20,
three instances of Lemma 12F22, and the cubic family of Lemma 12F15.
The stricter root-merger certificate of Lemma 12F25 certified 75,771
families.  Its deficit of 201 is exactly the number, through this cutoff, of
the 171 circuit-family instances of Lemma 12F18, the eighteen instances of
Lemma 12F24, and those same twelve one-root exceptions.
The theorem-backed `V_2` test, allowing either the original four roots or
one coalescence, certified 75,592 families.  Every one of its remaining 380
families belonged to Lemma 12F15, 12F18, 12F24A, or 12F24B; there were zero
unclassified cases.  Of the exact root-merger certificates, 197 lay outside
the `V_2` test, and all 197 were again in these symbolic families.
There were no classification failures.  Thus every unequal family through
label forty has either the jet certificate or one of the symbolic exceptional
proofs.  The exact chamber certificates above now replace this bounded audit
and discharge the formerly open integer comparison in Assumption 12F26A.

Two further support-only shortcuts were tested and rejected.  Prefix-cut
upper bounds predict three terms for source zero at `(3,5;2,2,2)`, whose
actual support has one.  Adding all interval-cut bounds and the deterministic
separated lower form predicts only one term for source one in the same
family, whose actual support has three.  Hence neither proposed support
polytope is exact.

When `q=r`, let `sigma` exchange the two equal odd terminals.  Both maps in
`(P12F.46)` commute with `sigma`, so `(P12F.47)` splits into its `sigma=+1`
and `sigma=-1` parts.  An antisymmetric invariant is divisible by `p_(qr)`:
it vanishes whenever the two terminal vectors are proportional, and
`p_(qr)` is the irreducible equation of that divisor.  After division the
quotient is `sigma`-symmetric and has terminal degree `q-1` at both
vertices.  The exact decomposition through label sixteen found no kernel in
the `sigma=+1` part and only

```text
dim ker(m)^-=dim ker(delta)^-=1
```

at `(P12F.48)`.  Hence the equal-terminal proof target is now: prove
injectivity of the symmetric part, and prove injectivity of the divided
antisymmetric part except for the displayed cubic allocation vector.

The symmetric terminal-tensor test has one uniform residual family, which
can be removed exactly.

**Lemma 12F27 (equal-terminal symmetric residual family).**  Let `q>=3` be
odd.  If

```text
q=r,                 (e_1,e_2,e_3)=(2,q+1,q+1),       (P12F.69)
```

then the `sigma=+1` merged source orbit sums are linearly independent.

**Proof.**  Write `a=e_1,b=e_2,c=e_3`, and over the fraction field of the
nonterminal bracket algebra write `A,B,C,V` for the four terminal linear
forms with roots `a,b,c,v`.  If `q>=5`, only `b` and `c` can be selected.
Substitution in `(P12F.50)`--`(P12F.53)` gives four orbit sums whose three
distinct symmetric terminal tensors are

```text
B^q odot (V C^(q-1)),       B^q odot C^q,
C^q odot (V B^(q-1)).                              (P12F.70)
```

These tensors are independent.  Indeed, put `B=X,C=Y,V=alpha X+beta Y`,
where `alpha beta!=0`.  The four forms

```text
B^q, C^q, V C^(q-1), V B^(q-1)
```

are independent: their coefficient supports contain respectively the four
distinguishing monomials `X^q,Y^q,X Y^(q-1),X^(q-1)Y`.  The tensors in
`(P12F.70)` are then three distinct basis monomials in the symmetric square
of their span.

The middle tensor occurs twice.  Its two nonterminal coefficients are

```text
p_(vb)p_(va)p_(ac),             p_(vc)p_(va)p_(ab),   (P12F.71)
```

which are not constant multiples.  Thus a ground-field relation first has
zero coefficients on the first and third tensors in `(P12F.70)`; its two
remaining constant coefficients vanish by `(P12F.71)`.

When `q=3`, selection of `a` adds two orbit sums.  Apart from the same
duplicated tensor, the five terminal tensors are

```text
(V A^2) odot (V B C),       (V A^2) odot (B C^2),
B^3 odot (V C^2),           B^3 odot C^3,
C^3 odot (V B^2).                                  (P12F.72)
```

Specialize `(V,A,B,C)` to the distinct affine slopes `(1,2,3,4)` and use
the coefficient basis
`X^3,X^2Y,XY^2,Y^3` in each tensor factor.  The minor of `(P12F.72)` in
flattened columns `0,1,2,3,5` has determinant `4392`.  Hence the five
generic tensors are independent.  The two copies of `B^3 odot C^3` are
again separated by `(P12F.71)`.  QED.

The strict checker `verify_su2_equal_terminal_exception.cpp` recomputes the
integer minor in the last paragraph.  A separate exact sweep through label
sixteen contained 343 nonempty equal-terminal families.  The full symmetric
terminal tensor had full rank in 336; the other seven were exactly
`(P12F.69)` and are covered by Lemma 12F27.  In the antisymmetric sector the
checker now uses the union of two sufficient tests: full antisymmetric
terminal-tensor rank, or full rank after division by `p_(qr)` and restriction
to `q=r`.  The division code explicitly discards terms with a remaining
`p_(qr)` factor before diagonal restriction and checks the resulting row
against the first-Wronskian identity.  This correction is essential: an
earlier diagnostic incorrectly deleted the remaining edge when merging the
terminals and therefore overstated the quotient rank.  Corrected exact spot
checks cover the cubic family and every previously identified residual.  A
clean current-source replay through label sixteen found 343 nonempty families
and certified all 343: the symmetric criterion passed all 343 after Lemma
12F27, and the corrected union criterion passed all 343 antisymmetric sectors.
The durable output is `certificates/su2_equal_terminal_corrected_16.log`; its
companion metadata records the source and binary hashes.  The remaining
equal-terminal task is the unbounded classification behind the corrected
union of sufficient tests.

The first half of that unbounded classification is now exact.

**Lemma 12F28A (exact equal-terminal strict-minimum classification).**  For
odd `q>=7`, failure of every strict single-root `V_2` filtration implies that
the ordered even triple is one of

```text
(2,q-1,q+1),       (2,q+1,q+1),
(q-1,q-1,q-1),     (q-1,q-1,q+1),
(q-1,q+1,q+1),     (q+1,q+1,q+1),
(q-1,q-1,2q-4),    (q-1,q-1,2q-2),
(q-1,q+1,2q-2),    (q-1,q+1,2q),
(q+1,q+1,2q-2),    (q+1,q+1,2q),
(q+1,q+1,2q+2).                                  (P12F.75)
```

For `q=3` the residual list has seven triples,

```text
(2,2,2), (2,2,4), (2,4,4), (2,4,6),
(4,4,4), (4,4,6), (4,4,8),
```

and for `q=5` it has twelve triples,

```text
(2,4,6), (2,6,6),
(4,4,4), (4,4,6), (4,4,8),
(4,6,6), (4,6,8), (4,6,10),
(6,6,6), (6,6,8), (6,6,10), (6,6,12).
```

**Certificate.**  For each of the four roots, the C++ verifier negates the
strict inequality of the two orbit-factor orders or the full `V_2` plus Hall
certificate of every leading-form group.  It then asserts that all four
roots fail and that the labels lie outside the displayed lists.  Z3's
quantifier-free linear-integer solver returns `unsat` over the entire
unbounded domain.  The durable transcript is
`certificates/su2_equal_strict_min_classification.log`; its metadata records
the clean source and executable hashes.  A separate concrete sweep through
label forty found exactly the same 174 residual instances and no others.

This lemma proves both terminal tensor sectors outside the displayed lists by
the ordinary one-root associated-graded argument.  It does *not* prove the
listed residuals.  In particular, valuations at distinct roots cannot be
lexicographically composed after passing to the first associated graded; an
experimental shortcut of that kind was rejected because it falsely predicted
full antisymmetric terminal rank in the cubic case.  The remaining valid task
is therefore an explicit tensor-minor proof for the families in `(P12F.75)`
and the nineteen fixed low-degree cases, with `(2,q+1,q+1)` already handled
by Lemmas 12F27--12F28 except for its fixed `q=3` antisymmetric minor.

The same residual family is also harmless in the divided antisymmetric
sector away from the cubic degree.

**Lemma 12F28 (antisymmetric residual family for `q>=5`).**  Let `q>=5` be
odd and assume `(P12F.69)`.  Divide each `sigma=-1` orbit difference by
`p_(qr)` and then set `q=r`.  The resulting four rows are linearly
independent.

**Proof.**  For degree-`q` binary forms put

```text
D(F,H)=F_X H_Y-F_Y H_X.
```

Up to a common nonzero scalar, diagonal restriction after division by
`p_(qr)` is `D(F,H)`.  The three distinct terminal Wronskians in this family
are

```text
D(B^q,V C^(q-1)),       D(B^q,C^q),
D(C^q,V B^(q-1)).                              (P12F.73)
```

Put `B=X,C=Y,V=alpha X+beta Y`, with `alpha beta!=0`.  The three expressions
in `(P12F.73)` are respectively

```text
q alpha(q-1) X^q Y^(q-2) + q^2 beta X^(q-1)Y^(q-1),
q^2 X^(q-1)Y^(q-1),
-q^2 alpha X^(q-1)Y^(q-1) - q beta(q-1)X^(q-2)Y^q.   (P12F.74)
```

Their outer monomials make them independent.  The middle Wronskian occurs
twice with opposite orientation; its two nonterminal coefficients are the
nonproportional monomials `(P12F.71)`.  The same ground-field coefficient
argument as in Lemma 12F27 separates those two rows.  QED.

The other twelve affine families in `(P12F.75)` admit a uniform exact minor.

**Lemma 12F29 (uniform tropical terminal minors).**  Let `q=2h+1>=7`.
For every family in `(P12F.75)` except `(2,q+1,q+1)`, both its symmetric
terminal tensors and its antisymmetric terminal tensors are linearly
independent.

**Proof.**  Write a terminal product form as

```text
M_a=V^(a_0) A^(a_1) B^(a_2) C^(a_3),
```

and, after a recorded permutation of the four generic roots, specialize

```text
(V,A,B,C)=(X,Y,X+Y,X+tY).                         (P12F.76)
```

For `0<=d<=q`, the coefficient of `X^(q-d)Y^d` in `M_a` is nonzero exactly
when

```text
a_1 <= d <= a_1+a_2+a_3,
```

and its highest `t`-degree is

```text
w_a(d)=min(a_3,d-a_1).                            (P12F.77)
```

Consequently the `t`-weight of a symmetric tensor entry in column `(d,e)`
is the maximum of

```text
w_a(d)+w_b(e),       w_b(d)+w_a(e).                (P12F.78)
```

The same formula applies to an antisymmetric entry whenever the two displayed
weights differ; then the leading terms cannot cancel.

The strict C++ checker chooses columns whose degrees are affine expressions
among constants, `q-c`, and `h+c`.  For each of the twelve families and each
sign it exhibits one determinant permutation having strictly greater total
`t`-weight than every competing permutation.  The checker does not infer
this from samples: after finding a candidate it constructs all entry weights
symbolically, audits the active source pattern, enumerates every determinant
permutation, and asks Z3's quantifier-free linear-integer solver whether the
chosen entries can disappear or a competitor can tie them for some `h>=3`.
All 24 negations are unsatisfiable.  Thus every displayed determinant is a
nonzero polynomial in `t`, so its generic determinant is nonzero.

The durable transcript is `certificates/su2_equal_tropical_uniform.log`.
Its metadata records the strict-build source, binary, and log hashes.  The
source-formula audit independently matched the symbolic exponent rows, with
multiplicity, against the noncrossing graph generator in all equal and
unequal tuples through label twelve.  QED.

**Corollary 12F30 (equal-terminal branch).**  The equal-terminal part of
`(P12F.47)` holds for every odd `q>=3` and every ordered positive even triple.

**Proof.**  Lemma 12F28A gives a strict one-root `V_2` filtration outside its
finite/affine residual lists.  For `q>=7`, Lemma 12F29 treats twelve affine
families and Lemmas 12F27--12F28 treat `(2,q+1,q+1)`.  For `q=3,5`, the exact
classification leaves only the nineteen listed triples.  The corrected
integer/modular rank certificate
`certificates/su2_equal_low_residuals.log` checks each of them: a full
terminal tensor or a full divided-diagonal quotient handles every sector
except the cubic antisymmetric sector, whose one-dimensional kernel is the
explicit cubic allocation relation of Lemma 12F15.  The special triple
`(3;2,4,4)` uses Lemma 12F27 in the symmetric sector and its certified full
divided-diagonal rank in the antisymmetric sector.  These are precisely the
kernel dimensions required in `(P12F.47)`.  QED.

### Arbitrary-background fixed-split diagnostic

The seven-factor identity `(P12F.47)` is not a formally multiplicative
statement.  The strict C++ mode
`--tau-even-background q r e_1 ... e_m` generalizes the merged `tau`-fixed
source to any `m>=3`.  It also constructs the complete negative
three-element-block source, the marked-`f` allocation map `iota_-`, the
global Plucker multiplication map `mu_-`, and the split-residual contraction
`C D`.  The latter is tested at the eight Vandermonde weights

```text
z_alpha=(mask(alpha)+1)^j,       1<=j<=8,
```

and at eight independent deterministic mask-weight vectors.  All maps are
constructed from noncrossing standard monomials; the displayed ranks are
over the certificate prime `1,000,000,007`.

Already with four even factors the stronger kernel equality fails while the
full allocation/payment dimensions still pass:

```text
[3,3,1,1,2,4,4,4]:
  dim W=52, ker(m)=6, ker(delta)=4;
  dim O_-=104, kappa_alloc=10, kappa_-=26;
  min_j dim ker(mu_- direct-sum C D_j)=6.
```

Thus arbitrary-background carrier relations occur even though all
seven-factor local families satisfy `(P12F.47)`.  More importantly, the
first fixed-split obstruction to the one-channel payment occurs at five
even factors:

```text
[3,3,1,1,4,4,4,4,4]:
  dim W=280, ker(m)=21, ker(delta)=6;
  dim O_-=560, kappa_alloc=20, kappa_-=121,
  kappa_car=101, dim Inv(V_3 tensor V_3 tensor V_4^tensor5)=196;
  dim ker(mu_- direct-sum C D_j)=29,       1<=j<=8.
```

All eight independent deterministic mask weights also give stacked-kernel
dimension `29`.  Hence the fixed `C D` family has a nine-dimensional
shortfall relative to `kappa_alloc` over the certificate field, despite the
large positive singlet reservoir.  This is not a counterexample to `Q3`:
for more than seven factors the signed subset expansion has positive supply
at several split sizes, whereas this diagnostic isolates only the
three-element-block source and one residual channel.  It instead rules out
the proposed proof that multiplies each seven-factor payment by an arbitrary
background monomial and treats every `3|(n-3)` layer independently.

The shortfall is presently a finite-field obstruction, not yet a rational
rank certificate: modular rank deficiency can in principle be caused by the
chosen prime.  An exact rational upper-rank certificate would make the
failure unconditional.  For roadmap purposes the correct next target is
therefore the cumulative prefix-Hall/outer-Turan formulation, or an enlarged
Plucker complex whose next chain group includes positive channels from all
split sizes.  A bounded-support local lift is no longer a viable default.

One stronger shortcut is false.  It is not always true that

```text
ker(mu_- direct_sum C D) subset ker(iota_-).
```

For `[3,3,1,1,2,2,2]`, the image under `iota_-` of the stacked kernel has
rank one.  Nevertheless `C D` has rank six on the seven-dimensional
negative kernel, while `kappa_car=5`, so `(P12F.35)` still holds.  The proof
must compare ranks; it cannot be reduced to this kernel inclusion.

The same weighted rank survives translated large-label families.  Direct
strict runs give

```text
[9,9,1,1,2,8,8]:    kappa_car=2, weighted rank=2, reservoir=25;
[11,11,1,1,2,10,10]: kappa_car=2, weighted rank=2, reservoir=31;
[9,11,1,1,2,10,10]:  kappa_car=2, weighted rank=2, reservoir=29.
```

Thus the observed mechanism is not confined to the low-label cutoff, but
these remain individual finite instances.

The ordinary `{a,f}` singlet contraction itself cannot be that map.  In
bracket coordinates it satisfies

```text
Omega_(af)(p_(af))=2,
Omega_(af)(p_(ai)p_(fj))=p_(ij),
```

but it acts on the resulting global invariant and therefore factors through
`mu_-`; it vanishes identically on `ker(mu_-)`.  The required payment map
must remember the factorization split before multiplication.

Thus the ordered edge-packet obstruction is already absent on the full
even-label, or `SO(3)`, subcone.  The interior `(OSD)` and higher-minus
steps remain separate requirements for a full `SO(3)` theorem.

A simpler checkerboard sign cone for the individual coefficients `f_Q(r,s)`
is false.  The exact diagnostic `search_su2_sp4_sign_cone.cpp` finds

```text
Q=[2,2,2,3,3],       f_Q(3,3)=37,
```

although the second index is odd.  Thus `(P12C.4)` is genuinely cumulative;
it cannot be replaced by the termwise rule `(-1)^s f_Q(r,s)>=0`.

The same two-packet split does *not* survive the affine wall term by term.
The finite-level checker `search_su2_level_edge_packet.cpp` first fails at
level `5`, word `[2,2,2,2,2,2,3,3]`, and wall target `5`: the `Y_2` packet
is `-17`, while the outer packet is `57`, leaving the positive full edge
coefficient `40`.  This does not contradict the successful finite
`(EDGE-CUT)` tests.  It proves that a uniform finite-level argument must
transfer reflected-wall supply from the outer packet, whereas the stable
ordinary proof may treat the two packets separately.

Ordinary total positivity of the individual factors is not available.  For
example, the `(2,0)|(2,0)` minor of `I+zN_2` is

```text
1+z-z^2.
```

Nor does the edge recurrence split into two separately nonnegative pieces.
Writing `E=E_P` and adjoining `p`, its new edge entry is

```text
E_(pP)(a,1)
 = sum_(r in a star p) E(r,1)
   +E(a,p-1)+E(a,p+1)-E(a-1,p)-E(a+1,p).
```

The last four-term expression is `-3` for prefix `[2,2,2]`, next label `3`,
and target `a=4`.  A related cut decomposition into a convolution of old
edge defects and the commutator `[N_1,C_P]m_suffix` also cannot be proved
termwise: the commutator is `-2` for word `[2,2,2,2,3]`, cut after the first
two factors, and target `4`.  The strict C++ modes
`search_su2_defect_cone.cpp --edge-step` and
`search_su2_edge_commutator.cpp` reproduce these exact obstructions.  Hence
the proof must retain the full exterior carrier through a crossing; local
Monge positivity loses necessary compensation.

A direct scalar four-functions shortcut does not prove the cut inequality.
The tempting local hypothesis

```text
m_S(a)m_T(b)
 <= m_(S intersect T)(0)
      sum_(r in a star b) m_(S union T)(r)
```

would imply boundary domination by the Reuter--Lovasz--Saks complementary-
subset inequality.  It is false even with disjoint target support.  Take the
indexed plus word `[1,1,2]`, let `S` contain the first `1` and the `2`, let
`T` contain the second `1` and the same `2`, and take `a=b=3`.  Both factors
on the left are one, while `S intersect T` is the singleton label `2` and
has invariant multiplicity zero.  Any four-functions argument must therefore
use a genuinely aggregated or representation-valued local inequality; scalar
log-supermodularity of the individual target multiplicities is insufficient.

Nor can the boundary reservoir in `(PH_t)` be replaced by the
never-updated term alone.  The `--exterior` mode of the same C++ diagnostic
tests this stronger assertion.  For six fundamental plus labels and target
`(2,2)`, the prefix through the fifth last-update class has saturated
residual `30`, whereas

```text
m_(1,1,1,1,1,1)(4)+m_(1,1,1,1,1,1)(2)=14.
```

A direct crystal-commutor shortcut also fails.  Group the factors assigned to
the two coordinates, place a highest element of weight `a` beside the lowest
element in a weight-`b` component, and use the type-`A_1` combinatorial
`R`-matrix to restore the original factor order.  The resulting element does
lie in a global component allowed by `a star b`, but the map loses the colour
assignment.  For plus word `[1,1,1,1,1,3]` and target `(2,2)`, the exact C++
diagnostic `character_ring_iter/search_su2_crystal_shuffle.cpp` sends 90
two-colour highest-weight sources to only four distinct images.  Hence a
crystal proof must switch the two colour carriers while retaining their
separate boundary endpoints; ordinary unshuffle-and-forget is insufficient.

Thus forgetting the two colours of the source path loses too much
information.  A successful injection must use boundary paths with a
nonempty invariant colour class; this is already forced in the
character-positive all-fundamental sector.

## A rank-one torus-square reduction for every even-minus sector

There is a second exact reduction which treats all even numbers of minus
factors at once and leaves only three Fourier charges.  For `j>0`, put

```text
C_j(u)=u^j+u^(-j),       S_j(u)=u^j-u^(-j),
W_p={p,p-2,...} intersect Z_(>0).
```

For each plus factor of label `p`, choose either `C_j` with `j in W_p`, or,
when `p` is even, choose the constant feature `1` with weight `2`.  For each
minus factor of label `p`, choose `S_j` with `j in W_p`.  A complete feature
choice `sigma` gives

```text
A_sigma(u)=product_(plus nonzero choices) C_j(u)
            product_(minus choices) S_j(u),
w_sigma=2^(number of zero plus choices).
```

When the number of minus factors is even, `A_sigma(u^(-1))=A_sigma(u)`.  Write

```text
a_r(sigma)=[u^r] A_sigma(u),
q(sigma)=a_0(sigma)^2-2a_2(sigma)^2
                         +a_0(sigma)a_4(sigma).
```

**Proposition 13 (three-charge torus-square identity).** For every signed
`SU(2)` character list with an even number of minus factors,

```text
J(M,P)=sum_sigma w_sigma q(sigma).                (P13.1)
```

**Proof.** Use torus variables `z,w` and the Weyl constant-term formula

```text
integral_SU(2) f = CT_z (1-z^2)f(z).
```

Since every character expression is invariant under inversion of either
torus variable, `1-z^2` may equivalently be replaced by

```text
omega(z)=1-(z^2+z^(-2))/2=-(z-z^(-1))^2/2.
```

Make the injective exponent change `z=uv`, `w=u/v`; a Laurent monomial has
zero `z,w` exponent if and only if its image has zero `u,v` exponent.  Direct
expansion of the `SU(2)` weight strings gives

```text
chi_p(uv)+chi_p(u/v)
 = sum_(j in W_p) C_j(u)C_j(v) + 2 indicator_(p even),

chi_p(uv)-chi_p(u/v)
 = sum_(j in W_p) S_j(u)S_j(v).                  (P13.2)
```

Thus the whole signed product is the positive feature sum

```text
sum_sigma w_sigma A_sigma(u)A_sigma(v).
```

The two Weyl factors become

```text
omega(uv)omega(u/v)
 = (1/4)[C_2(u)-C_2(v)]^2.                       (P13.3)
```

For even minus parity, inversion symmetry gives `a_(-r)=a_r`.  Taking the
constant term of `(P13.3) A_sigma(u)A_sigma(v)` now gives exactly

```text
a_0^2-2a_2^2+a_0 a_4.
```

Summing the features proves `(P13.1)`.  QED.

The feature sum has an intrinsic Hilbert-space form.  For every label `p`,
let `H_p` have the orthonormal basis

```text
e_m,       m=p,p-2,...,-p,
```

and let `R_p e_m=e_(-m)`.  Put `P_p^+=(I+R_p)/2`,
`P_p^-=(I-R_p)/2`, and

```text
w_p(u)=sum_m u^m e_m,
v_p^epsilon(u)=sqrt(2) P_p^epsilon w_p(u).
```

For a signed factor list, form the tensor-valued Laurent polynomial

```text
V(u)=tensor_i v_(p_i)^(epsilon_i)(u)=sum_r u^r v_r.
```

**Proposition 14 (three-charge Gram form).**  In the notation of
Proposition 13,

```text
<v_r,v_s> = sum_sigma w_sigma a_r(sigma)a_s(sigma)=B_(r,s),

J(M,P)=||v_0||^2-2||v_2||^2+<v_0,v_4>.          (P14.1)
```

**Proof.**  For every positive weight `j`, use the normalized reflection
eigenvectors

```text
f_j^+=(e_j+e_(-j))/sqrt(2),
f_j^-=(e_j-e_(-j))/sqrt(2).
```

Then

```text
v_p^+(u)=sum_(j in W_p) C_j(u)f_j^+
          +sqrt(2) indicator_(p even)e_0,

v_p^-(u)=sum_(j in W_p) S_j(u)f_j^- .            (P14.2)
```

The tensor products of the displayed local vectors form an orthonormal
feature basis.  The component indexed by `sigma` in `V(u)` is exactly
`sqrt(w_sigma) A_sigma(u)`: each zero plus choice contributes the factor
`sqrt(2)`, and every other choice has unit norm.  Comparing the coefficient
of `u^r` proves the Gram identity.  Substitution in `(P13.1)` gives
`(P14.1)`.  QED.

Thus the all-even-minus problem is equivalently the single structured norm
inequality

```text
2||v_2||^2 <= ||v_0||^2+<v_0,v_4>.              (P14.3)
```

The qualification "structured" is essential: a general positive
semidefinite Gram matrix does not satisfy `(P14.3)`.  Here the vectors are
the charge projections of tensor products of reflection eigenspaces of
complete `SU(2)` weight strings.  A successful operator proof must use that
complete-string and nested-weight structure, rather than positivity of the
Gram matrix alone.

There is a further exact telescoping which moves the three-charge expression
to the outer edge.  Remove one signed factor of label `p`, write its sign as
`epsilon=+1` for a plus factor and `epsilon=-1` for a minus factor, and call
the remaining signed product `R(x,y)`.  Because the original number of minus
factors is even,

```text
R(y,x)=epsilon R(x,y).                              (P15.1)
```

Let

```text
G_R(x)=integral R(x,y) dmu(y)=sum_(r>=0) g_r chi_r(x).
```

For the feature choices `tau` of the remaining list, write

```text
B_tau(u)=sum_r b_r(tau)u^r,
b_(-r)(tau)=epsilon b_r(tau),

H_r(R)=sum_tau w_tau[
          b_r(tau)^2-b_(r-2)(tau)b_(r+2)(tau)].         (P15.2)
```

When `epsilon=-1`, the convention in `(P15.2)` in particular gives
`b_0=0`.

**Proposition 15 (outer Turan reduction).**  With the preceding notation,

```text
J(M,P)=2g_p=2[H_p(R)-H_(p+2)(R)].               (P15.3)
```

Consequently the full central-character `Q3` theorem is equivalent to the
following outer-tail statement:

```text
If every label occurring in R is at most p and
R(y,x)=epsilon R(x,y), then

H_p(R)>=H_(p+2)(R).                              (OTM)
```

Indeed, for any nonempty signed list with an even number of minus factors,
remove a factor having a largest label.  Its sign is exactly the exchange
parity of the remainder, so `(OTM)` and `(P15.3)` apply.  Conversely, adding
a factor of label `p` and sign `epsilon` to any remainder in `(OTM)` gives
an even-minus signed list, so `Q3` implies the displayed inequality.

**Proof.**  By `(P15.1)` and exchange invariance of product Haar measure,

```text
integral integral [chi_p(x)+epsilon chi_p(y)]R(x,y)
 =2 integral chi_p(x)G_R(x)dmu(x)=2g_p.           (P15.4)
```

It remains to identify the coefficient.  Fix one remaining feature
polynomial `B(u)=sum_r b_r u^r`.  For a nonzero local feature `j` of the
removed factor, multiply `B` by `C_j` when `epsilon=+1` and by `S_j` when
`epsilon=-1`, and insert its coefficients at charges `0,2,4` into the
quadratic form in Proposition 13.  Sum over

```text
j=p,p-2,...,1 or 2,
```

including the weight-two zero feature when `epsilon=+1` and `p` is even.
All interior terms cancel, leaving exactly

```text
2[(b_p^2-b_(p-2)b_(p+2))
   -(b_(p+2)^2-b_p b_(p+4))].                   (P15.5)
```

For `p=1`, the first lower coefficient is
`b_(-1)=epsilon b_1`; for even `p` and `epsilon=-1`, the apparent zero-charge
boundary terms vanish because `b_0=0`.  Thus `(P15.5)` is valid uniformly.
Summing it with the remaining feature weights gives the second equality in
`(P15.3)`.  Equivalently, applying the Weyl constant-term operator first in
`y` shows that `H_r(R)` is the torus weight coefficient of charge `r` in
`G_R`; the usual difference of adjacent weight coefficients gives
`g_r=H_r(R)-H_(r+2)(R)`.  QED.

The strict symbolic C++ diagnostic
`character_ring_iter/inspect_su2_torus_boundary.cpp` expands the local
quadratic form and displays the cancellation in `(P15.5)` for either sign.
The new target `(OTM)` is still unproved, but it has two advantages over the
raw three-charge inequality: its charge lies at or beyond every individual
remaining label, and it asks for monotonicity only at that outer edge, not
at every charge.

Outer Turan monotonicity is not featurewise.  The strict C++ search
`character_ring_iter/search_su2_turan_atom.cpp` finds, at outer label `3`,

```text
B=S_1^3 C_1^4 C_2,       H_3(B)=4,       H_5(B)=6.
```

Even restricting every magnitude to be even does not repair the atomic
claim: at outer label `2`, `B=C_2^2` gives `H_2=-2` and `H_4=1`.  These
failures do not contradict `(OTM)`, because a complete irreducible weight
string sums over all its nested features; for example an even plus factor
also has its weight-two zero feature.  They show that a proof must preserve
those complete-string groupings.  The symbolic diagnostic
`character_ring_iter/inspect_su2_turan_update.cpp` expands the corresponding
one-factor update and can be used to search for a closed cone of grouped
outer Turan forms.

The outer reduction proves a new unbounded near-extremal region.  Continue
to write `q_1,...,q_n` for the positive labels in the remainder, put
`S=sum_i q_i`, and let `s_i` be `+1` or `-1` according as its factor is plus
or minus.

**Proposition 16 (outer two-layer theorem).**  If the removed largest label
`p` satisfies

```text
p>=S-2,
```

then its signed-product integral is nonnegative, with no bound on `n` or on
the labels.  More precisely, the outer coefficient `g_p` in Proposition 15
is zero when `p>S` or `p` has parity opposite to `S`, is one when `p=S`, and
when `p=S-2` it is

```text
g_(S-2)=n-1+sum_(i<j, q_i=q_j=1) s_i s_j >=0.    (P16.1)
```

**Proof.**  Expand the remainder by choosing a subset `T` of its factors for
the `y` coordinate.  If `N(T)` is the invariant multiplicity of those
labels, then

```text
g_p=sum_T (product_(i in T)s_i) N(T)
             multiplicity_(V_p)(tensor_(i notin T)V_(q_i)).       (P16.2)
```

The tensor product in the last factor has highest weight `S-sum_(i in T)q_i`.
This immediately gives the assertions for `p>S`, opposite parity, and
`p=S`.  For `p=S-2`, the empty subset contributes the multiplicity of the
next-to-highest constituent in a product of `n` positive-label irreducibles.
That multiplicity is `n-1`: the weight of height two below the top has
multiplicity `n`, while the unique top constituent accounts for one copy.

A nonempty invariant tensor has even total label at least two.  At total
label two it consists of exactly two fundamental labels, its invariant
multiplicity is one, and the complementary tensor product reaches `V_(S-2)`
through its unique top constituent.  Any subset of larger total label leaves
complementary highest weight below `S-2`.  This proves the equality in
`(P16.1)`.

If `r` of the `q_i` are fundamental and their signs have sum `d`, then

```text
sum_(i<j, q_i=q_j=1) s_i s_j=(d^2-r)/2>=-floor(r/2).
```

Since `r<=n`, one has `n-1-floor(r/2)>=0` (including `n=1`), proving the
claim.  QED.

The next layer can also be evaluated and bounded uniformly.

**Proposition 17 (outer four-defect theorem).**  In the notation of
Proposition 16, the outer coefficient is nonnegative also when

```text
p=S-4>=max_i q_i.
```

Hence every signed-product integral is nonnegative whenever a largest
selected label satisfies `p>=S-4`.

**Proof.**  Put `r=#{i:q_i=1}`, `t=#{i:q_i=2}`, and let `h=n-r-t` count the
labels at least three.  Write

```text
e_k=[z^k] product_(q_i=1)(1+s_i z),
y=sum_(q_i=2)s_i,
f_2=sum_(i<j, q_i=q_j=2)s_i s_j.
```

For the empty `y`-subset, the multiplicity of `V_(S-4)` is

```text
A=n(n-3)/2+(t+h).                                  (P17.1)
```

Indeed, at weight height four below the top one may lower two distinct
factors once, or lower once more any factor of label at least two.  The
weight multiplicity is therefore `binomial(n,2)+(t+h)`; subtracting the top
constituent and the `n-1` next-to-top constituents gives `(P17.1)`.

Every nonempty invariant subset that can contribute has total label two or
four.  The possibilities and their invariant multiplicities are

```text
[1,1]: 1,       [2,2]: 1,       [1,1,2]: 1,
[1,1,1,1]: 2.
```

The label-two subsets leave a next-to-top multiplicity `n-3`; the label-four
subsets leave the unique top constituent.  Formula `(P16.2)` consequently
becomes

```text
g_(S-4)=A+(n-3+y)e_2+f_2+2e_4.                  (P17.2)
```

It remains to prove that this signed expression is nonnegative.  First fix
the signs of the `r` fundamental factors.  Adding a plus label `2` changes
the right side of `(P17.2)` by

```text
n+2e_2+y=d^2+h+2t_+ >=0,
```

where `d` is the sum of the fundamental signs and `t_+` is the number of
plus label-two factors already present.  Adding a minus label `2` changes it
by

```text
n-y=r+h+2t_- >=0.
```

Adding a label at least three changes it by

```text
n+e_2=(r+2t+2h+d^2)/2 >=0.
```

Thus `(P17.2)` is minimized, for fixed fundamental signs, when every
remaining label is fundamental.  In that case `n=r`, and Newton's identities
give

```text
e_2=(d^2-r)/2,
e_4=[d^4+(8-6r)d^2+3r^2-6r]/24.
```

Substitution in `(P17.2)` cancels the other terms and leaves

```text
g_(r-4)=[d^4-10d^2+3r^2-6r]/12.                (P17.3)
```

Admissibility of `p=r-4>=1` gives `r>=5`.  For `r=5`, the allowed values
`|d|=1,3,5` give respectively `3,3,35`.  For `r>=6`, the term
`d^4-10d^2` is at least `-25`, while `3r^2-6r>=72`; hence `(P17.3)` is
positive.  This proves the claim.  QED.

The independent strict C++ verifier
`character_ring_iter/verify_su2_outer_layers.cpp` constructs the signed
two-coordinate fusion product, extracts `g_(S-4)`, and compares it with
`(P17.2)`.  It passed all 60,972 admissible signed remainders with labels at
most five and at most seven factors.

Thus any still-open configuration must lie at least six weight units inside
the outer tensor-product boundary: for a largest selected label `p`, it must
have `p<=S-6` with the matching parity.

The proof of Proposition 17 does not extend by simple monotonicity under
adjoining complete factors.  The strict exact diagnostic
`character_ring_iter/search_su2_outer_append.cpp` first fails at deficit six:
five plus label-two factors have outer coefficient `45`, while adjoining one
minus label-two factor changes the corresponding fixed-deficit coefficient
to `44`.  Both values are positive; the failure only rules out an induction
which demands that every appended factor increase the coefficient.

The deficit-six layer nevertheless has a finite exact formula.  For labels
`1`, `2`, and `3`, let `e_j^(a)` be the degree-`j` elementary symmetric sum
of the signs of the factors having label `a`.  Put

```text
r=#{i:q_i=1},       t=#{i:q_i=2},       n=#factors,

A=binomial(n+1,3)-r(n-1)-t,
B=(n-2)(n-5)/2+n-r.
```

**Lemma 17A (exact six-defect formula).**  If `p=S-6>=max_i q_i`, then

```text
g_(S-6)
 = A+B e_2^(1)
   +(n-3)e_2^(2)+(n-4)e_2^(1)e_1^(2)
   +2(n-5)e_4^(1)+e_2^(3)
   +e_1^(1)e_1^(2)e_1^(3)+e_3^(2)
   +e_3^(1)e_1^(3)+2e_2^(1)e_2^(2)
   +3e_4^(1)e_1^(2)+5e_6^(1).                  (P17A.1)
```

**Proof.**  In `(P16.2)`, a contributing contracted subset has total label
at most six.  The empty subset leaves the third constituent below the top.
The weight coefficient at lowering degree three is

```text
binomial(n+2,3)-rn-t,
```

while the degree-two weight coefficient is
`binomial(n+1,2)-r`; their difference is `A`.

The total-label-two subset `[1,1]` leaves a deficit-four complementary
coefficient `B`.  At total label four, the invariant subsets and
multiplicities are

```text
[2,2]:1,       [1,1,2]:1,       [1,1,1,1]:2,
```

and their complements contribute at the next-to-top layer.  At total label
six, the complete list is

```text
[3,3]:1,             [1,2,3]:1,       [2,2,2]:1,
[1,1,1,3]:1,         [1,1,2,2]:2,
[1,1,1,1,2]:3,       [1,1,1,1,1,1]:5.
```

Their complements contribute through the unique top constituent.  Summing
the signed choices in these lists gives exactly `(P17A.1)`.  QED.

After Proposition 5, every occurrence of a fixed label has one sign.  Write
those three signs as `alpha,beta,gamma`.  Then

```text
e_j^(1)=alpha^j binomial(r,j),
e_j^(2)=beta^j binomial(t,j),
```

and similarly for label three.  The sign-dependent part of `(P17A.1)` uses
only

```text
beta,       delta=alpha gamma.
```

Thus only four sign patterns remain.  The case `(beta,delta)=(+1,+1)` is
the all-plus/central-parity-aligned case already covered by Proposition 20.

**Proposition 17B (six-defect positivity).**  Under the hypotheses of Lemma
17A,

```text
g_(S-6) >= 0.                                             (P17B.1)
```

**Proof.**  First remove all factors of label at least four.  Re-adjoining
one such factor changes the six-defect coefficient by a sum of the outer
coefficients at deficits zero, two, and four.  Those coefficients are
nonnegative by Propositions 16 and 17.  It is therefore enough to treat
labels one, two, and three.

Write

```text
R_j=binomial(r,j),   T_j=binomial(t,j),   U_j=binomial(u,j),
C(r)=-R_2-2R_3+4R_4+10R_5+5R_6,
```

where `u` is the number of label-three factors.  Substitution in `(P17A.1)`
is most transparent in the Newton basis.  If `(beta,delta)=(+1,-1)`, it
gives

```text
g_(S-6)=H_+(r,t)+K_+(r,t,u),

H_+=5T_3+(2r+5R_2)T_2+(-1+6R_3+5R_4)t+C(r),              (P17B.2)
```

where

```text
K_+=2U_2+tu+tU_2+2T_2u+rU_2
    +R_2u+R_2U_2+2R_2tu+2R_3u+2R_4u >= 0.                (P17B.3)
```

For `r>=5`, the identity

```text
C(5+x)=70 binomial(x,1)+154 binomial(x,2)
       +168 binomial(x,3)+104 binomial(x,4)
       +35 binomial(x,5)+5 binomial(x,6)                  (P17B.4)
```

shows that every term of `H_+` is nonnegative.  For `0<=r<=4`, the explicit
cubics are

```text
r=0: 5T_3-t,
r=1: 5T_3+2T_2-t,
r=2: 5T_3+9T_2-t-1,
r=3: 5T_3+21T_2+5t-5,
r=4: 5T_3+38T_2+28t-10.
```

They are negative only at

```text
(r,t)=(0,1),(0,2),(1,1),(2,0),(2,1),(3,0),(4,0).          (P17B.5)
```

Admissibility forces respectively

```text
u>=3,2,2,3,2,2,2.
```

At those lower endpoints, the triples
`(g,Delta_u g,Delta_u^2 g)` are respectively

```text
(12,13,6), (10,13,6), (5,10,6), (18,19,8),
(12,17,8), (13,22,10), (34,41,14).                        (P17B.6)
```

They are nonnegative, and the degree in `u` is two.  This proves the
`(+1,-1)` case.

It remains to take `beta=-1`; the two values of `delta` have the common
`u`-free part

```text
H_-(r,t)=3T_3+(R_2+2r)T_2+(2R_2-R_4-1)t+C(r).             (P17B.7)
```

For `delta=+1` the remainder is

```text
K_-^+=2U_2+tu+tU_2+2T_2u+rU_2
      +R_2u+R_2U_2+4R_3u+2R_4u,
```

and for `delta=-1` it is

```text
K_-^-=2U_2+tu+tU_2+2T_2u+rU_2+2rtu
      +R_2u+R_2U_2+2R_3u+2R_4u.                          (P17B.8)
```

Both are nonnegative.  For `0<=r<=7`, expand `(P17B.7)` in the Newton basis
in `t` at the bases

```text
t_0(r)=4,4,3,3,2,2,1,0.
```

The coefficients are

```text
r       Newton coefficients at t_0(r)
0       8,17,12,3
1       20,25,14,3
2       20,25,14,3
3       40,41,18,3
4       24,41,20,3
5       48,57,26,3
6       84,41,30,3
7       294,6,35,3
```

Their minimum is `3`, so these expansions cover every `t>=t_0(r)`.  Direct
evaluation below the eight bases gives a negative result only
at

```text
(r,t)=(0,1),(0,2),(1,1),(2,0),(3,0),(4,0).                (P17B.9)
```

Admissibility forces `u>=3,2,2,3,2,2`.  At these endpoints the triples
`(g,Delta_u g,Delta_u^2 g)` are

```text
delta=+1: (12,13,6), (10,13,6), (5,10,6),
          (18,19,8), (17,24,10), (50,49,14),

delta=-1: (12,13,6), (10,13,6), (9,12,6),
          (18,19,8), (13,22,10), (34,41,14).              (P17B.10)
```

Again all later `u` are covered by the Newton expansion.

Finally suppose `r>=8`.  Dropping the nonnegative term `3T_3` from
`(P17B.7)` leaves the quadratic

```text
Q(t)=a binomial(t,2)+ell t+C(r),
a=R_2+2r,                ell=2R_2-R_4-1.
```

Its real minimum is nonnegative if

```text
D(r)=8aC(r)-(2ell-a)^2 >= 0.                              (P17B.11)
```

The exact Newton expansion at `r=8+x` has coefficients

```text
290204, 538260, 684027, 618414, 399894,
181560, 55140, 10080, 840,                                (P17B.12)
```

in degrees zero through eight.  Hence `D(r)>=0` for every integer `r>=8`.
This proves both `beta=-1` cases.  Together with the
`(beta,delta)=(+1,+1)` case, all four possibilities are covered.  QED.

The strict verifier `character_ring_iter/verify_su2_outer_layers.cpp`
compared `(P17A.1)` with the full signed fusion calculation in 60,080 cases
through labels five and seven factors.  The independent formula analyzer
`character_ring_iter/analyze_su2_outer_six_formula.cpp` checked 1,947,556
integer sign/count cases through thirty factors; its minimum was `3`.  It
also checks every finite Newton coefficient in `(P17B.4)`--`(P17B.12)` using
exact integer arithmetic and reports `SU2_OUTER_SIX_UNIFORM_CERTIFICATE
PASS`.

Thus any still-open configuration must lie at least eight weight units
inside the outer tensor-product boundary: for a largest selected label `p`,
it must have `p<=S-8` with matching parity.

The next layer has a finite exact reduction but is not yet proved in all
sign sectors.  At deficit eight, any factor of label at least five can be
removed: re-adjoining it changes `g_(S-8)` by the sum of the already proved
outer coefficients at deficits zero, two, four, and six.  After the
central-parity sign normalization, the remaining labels `1,2,3,4` give only
eight uniform sign patterns.  The admissible count domain is the upward
closure of 37 minimal bases.

The exact C++ analyzer
`character_ring_iter/analyze_su2_outer_eight.cpp` computes `g_(S-8)` from
subset invariant multiplicities using arbitrary-precision integers.  It
passed all 370,224 admissible count/sign cases through thirty factors, with
minimum `7`.  Its multivariate Newton certificate uses the exact degree
bound `(8,4,4,4)`.  Its total count-degree is at most eight: a contracted
subset of total label `2j<=8` contributes a binomial selection polynomial
of degree at most `2j`, while the complementary constituent lies
`4-j` lowering steps below its top and has count-degree at most `4-j`;
the sum is at most `4+j<=8`.  Consequently every transformed ray direction
also has degree at most eight.  The all-plus pattern and the normalized pattern

```text
(s_1,s_2,s_3,s_4)=(+1,+1,+1,-1)
```

have nonnegative Newton coefficients on all 37 minimal orthants, so those
two sectors are proved for unbounded counts.  The original single-orthant
certificate fails, not the positivity claim, in the other six sectors:
patterns two through five have coefficient `-1` at base `(0,0,0,3)` and
Newton degree `(1,1,2,0)`, while patterns six and seven first fail at degree
`(4,1,0,1)`.

Two of those six sectors admit a finite recursive orthant certificate.  If
a Newton expansion at a base `b` has a negative coefficient, choose one
free count coordinate `x` and partition the orthant into the slice `x=b_x`
and the shifted tail `x>=b_x+1`.  A nonnegative Newton expansion on every
leaf proves positivity on the original orthant by induction on this binary
partition.  The `recursive` mode performs this exact partition using
arbitrary-precision integers and the degree bound `(8,4,4,4)`.

**Proposition 17C (half of the eight-defect sign sectors).**  Under the
deficit-eight hypotheses, `g_(S-8)>=0` whenever the normalized sign on every
label-two factor is positive.  Equivalently, sign patterns zero through
three are nonnegative for all admissible counts.

**Proof.**  Patterns zero and one have the direct 37-orthant certificates
above.  For each of patterns two and three, the recursive verifier covers
the same 37 minimal orthants with 71 exact nodes, 54 nonnegative Newton
leaves, and 17 binary splits; the maximum recursion depth is four.  The
slice/tail partition is exhaustive at every split, and every leaf expansion
has nonnegative integer coefficients, proving the assertion on the upward
closure of all 37 bases.  This upward closure is precisely the admissible
count domain.  QED.

Pattern four, with normalized signs `(s_1,s_2,s_3,s_4)=(+,-,+,+)`, has a
ratio-sensitive certificate.  Write `(r,t,u,v)` for the numbers of factors
of labels one through four.  Split the `(r,t)` quadrant into the three
unimodular chambers

```text
lower:   (r,t)=(a+b,b),
middle:  (r,t)=(a+b,a+2b),
upper:   (r,t)=(n,2n+m).                         (P17D.1)
```

The admissible set is upward closed in the four coordinates of each
chamber.  Its product-order minimal sets have respectively 32, 20, and 20
elements.  This enumeration is exhaustive with transformed coordinates at
most 16.  Indeed every chamber direction has positive total-label weight at
most five.  If a minimal point had total label at least 17, removing any
present direction would leave total label at least 12, hence would remain
admissible because the largest available label is four, a contradiction.

All 32 lower-chamber minimal orthants have nonnegative Newton coefficients.
In the middle chamber, ten of the twenty minimal orthants are likewise
coefficientwise nonnegative.  Each exceptional polynomial coefficientwise
dominates either a shift of `H_0` or `H_1`, where, with `N=X+Y`,

```text
H_0(X,Y)=sum_(i=0)^6 (i^2-12i+52) C(X,i)C(Y,6-i)
        +sum_(i=0)^7 (4i-21)      C(X,i)C(Y,7-i)
        +14 sum_(i=0)^8           C(X,i)C(Y,8-i),        (P17D.2)

H_1(X,Y)=sum_(i=0)^6 (i^2-8i+26) C(X,i)C(Y,6-i)
        +sum_(i=0)^7 (4i-7)       C(X,i)C(Y,7-i)
        +14 sum_(i=0)^8           C(X,i)C(Y,8-i).        (P17D.3)
```

Vandermonde's identity and its first two factorial moments give, for
`N>=6`,

```text
H_j(X,Y)=C(N,6) P_j(N,X)/(4N(N-1)),               (P17D.4)

P_0=120X^2+8(2N^2-47N+30)X
             +N(N-1)(N^2-25N+322),

P_1=120X^2+8(2N^2-35N+18)X
             +N(N-1)(N^2-17N+170).
```

The discriminants of these quadratics in `X` are

```text
disc_X(P_0)=-16(14N^4-28N^3+1094N^2+1620N-3600),
disc_X(P_1)=-16(14N^4+20N^3+422N^2-60N-1296).     (P17D.5)
```

Both parenthesized polynomials are positive for `N>=6`; the leading
coefficient 120 is positive.  Hence `H_0,H_1>0` there, while both vanish
when all binomial degrees exceed `N`.  This proves every middle exception,
including the shifted `H_0` cases.

For the upper chamber, every minimal orthant has only the four possible
negative Newton degrees

```text
(4,1,0,1), (4,1,1,0), (6,1,0,0), (7,0,0,0).      (P17D.6)
```

The first two coefficients are `-1`, the third is `-4`, and the last is
`-21`, `-7`, or nonnegative.  Reserve coefficient two at degrees
`(4,2,0,0)`, `(4,0,0,2)`, `(4,0,1,1)`, and `(4,0,2,0)`, as well as three
units at `(4,0,0,0)`.  Put `x=n-4`, and denote the remaining pure and
linear-in-`m` coefficients by `c_(i,0)` and `c_(i,1)`.  After factoring
`C(n,4)`, define

```text
A(x)=sum_(j=0)^4
       (c_(4+j,0)-3 delta_(j,0))/C(4+j,j) C(x,j),
L(x)=sum_(j=0)^2 c_(4+j,1)/C(4+j,j) C(x,j).       (P17D.7)
```

The pure block is

```text
B=A(x)+L(x)m+(2/3)m(m-1).                         (P17D.8)
```

With `Ahat=210A` and `Phat=210L-140`,

```text
210B=140m^2+Phat(x)m+Ahat(x).                     (P17D.9)
```

The exact verifier checks, separately for all twenty upper minimal bases,
that `Ahat>=0` on every nonnegative integer `x`.  It finds the first point
`T` at which `Phat` is negative, proves `Phat` remains decreasing there,
and certifies

```text
560 Ahat(x)-Phat(x)^2>=0  for x>=T.               (P17D.10)
```

Here `T` ranges only from 30 through 35.  Finite prefixes are evaluated
exactly, and each infinite tail is proved by a nonnegative shifted Newton
expansion of degree at most four.  Thus `(P17D.9)` is nonnegative whether
`Phat` is nonnegative or negative.

The remaining mixed block, with `s=u+v`, is

```text
C(n,4){3+(1/3)m(m-1)+s(s-1)-ms}.                 (P17D.11)
```

After multiplication by three, its discriminant as a quadratic in `m` is
`-3(s-3)^2-8<0`.  It is therefore positive.  Every unused Newton
coefficient is nonnegative.

**Proposition 17D (a fifth eight-defect sign sector).**  Under the
deficit-eight hypotheses, `g_(S-8)>=0` for normalized sign pattern four,
`(s_1,s_2,s_3,s_4)=(+,-,+,+)`.

**Proof.**  The three chambers in `(P17D.1)` cover the nonnegative
`(r,t)` quadrant.  The exhaustive exact certificate proves all lower
orthants directly, all middle exceptions by `(P17D.2)--(P17D.5)`, and all
upper exceptions by `(P17D.6)--(P17D.11)`.  QED.

The certificate is replayed by

```text
./analyze_su2_outer_eight 1 pattern4-chambers 4
```

and is included in `certificates/su2_outer_eight_all_patterns.log`.  It uses
arbitrary-precision integers throughout and reports
`PATTERN4_GLOBAL_CHAMBER_CERTIFICATE PASS`.

Changing only the normalized label-four sign gives pattern five.  On every
transformed minimal orthant its polynomial has the same negative template
as pattern four; the coefficients affected by the sign change remain in the
nonnegative residual.  The exact dominance, middle discriminants, and all
twenty upper integer-tail certificates replay without alteration.

**Proposition 17E (the sixth eight-defect sign sector).**  Under the
deficit-eight hypotheses, `g_(S-8)>=0` for normalized sign pattern five,
`(s_1,s_2,s_3,s_4)=(+,-,+,-)`.

**Proof.**  Apply the exhaustive three-chamber certificate of Proposition
17D with the pattern-five coefficients.  All direct leaves pass; every
exceptional polynomial dominates the same nonnegative analytic template.
QED.

Patterns six and seven have negative label-three sign and add one new
negative family.  In the middle chamber at the only four affected minimal
bases, put `N=a+b` and

```text
R(a,b)=sum_(i=0)^5 (i-4)C(a,i)C(b,5-i)
      =a C(N-1,4)-4C(N,5).                         (P17F.1)
```

The full polynomial coefficientwise dominates the old `H_0` template plus

```text
J-H_0=C(N,4)+R(a,b)u+2C(N,4)C(u,2).               (P17F.2)
```

For `N<4` this block vanishes.  At `N=4` it is `1+u(u-1)`.  At `N=5` it is

```text
5u^2+(a-9)u+5;
```

its values at `u=0,1` are `5,a+1`, and it is increasing thereafter.  For
`N>=6`, use `P_0` from `(P17D.4)` and set

```text
Q=5a(N-4)-4N^2+11N,

F=5N(N-4)(N-5)P_0-6(N-1)Q^2+600N^2(N-1).         (P17F.3)
```

After division by `C(N,4)`, the discriminant margin of `(P17F.2)` together
with `H_0` is

```text
4(1+H_0/C(N,4))-(R/C(N,4)-1)^2
   =F/(150N^2(N-1)).                               (P17F.4)
```

The polynomial `F(a,b)` has total degree at most seven.  Its Newton
expansion is nonnegative on each of the seven orthants based at
`(i,6-i)`, `0<=i<=6`.  These orthants are exactly `a+b>=6`, so `(P17F.4)`
is nonnegative everywhere required.  This closes the new middle family.

For the upper chamber, write `x=n-4`.  At each of the twenty minimal bases,
let `c_(i,j,k)` be the exact Newton coefficient at degree `(i,j,k,0)` and
define

```text
A(x)=sum_(d=0)^4 c_(4+d,0,0)/C(4+d,d) C(x,d),
L(x)=sum_(d=0)^2 c_(4+d,1,0)/C(4+d,d) C(x,d),
G(x)=sum_(d=0)^2 c_(4+d,0,1)/C(4+d,d) C(x,d).      (P17F.5)
```

The selected template, after factoring `C(n,4)`, is

```text
A+m(m-1)+u(u-1)+v(v-1)+2uv-mu-mv+Lm+Gu.          (P17F.6)
```

Put `s=u+v`.  For fixed `s`, the minimum of `Gu` over `0<=u<=s` is
`min(G,0)s`.  Thus `(P17F.6)` is bounded below by

```text
A+m^2+s^2-ms+pm+qs,
p=L-1,  q=-1+min(G,0).                             (P17F.7)
```

Scale `Ahat=210A`, `phat=210p`, and `qhat=210q`.  The verifier proves that
`G>=0` permanently after `x=0` or `x=3`, depending on the minimal base.
On that tail `qhat=-210`, and the positive-definite quadratic in `(m,s)` is
nonnegative whenever

```text
630 Ahat-(phat^2+phat qhat+qhat^2)>=0.             (P17F.8)
```

For every base, the degree-four polynomial in `(P17F.8)` has a nonnegative
shifted Newton expansion from some `x<=13`.  The preceding finite values of
`x` are certified by exact integer minimization.  This minimization is
exhaustive: with `M=abs(phat)+abs(qhat)` and
`rho=sqrt(m^2+s^2)`,

```text
210(m^2+s^2-ms)+phat m+qhat s >= 105rho^2-M rho.
```

Consequently a minimizer lies in the explicitly enumerated square
`0<=m,s<=floor(M/105)+2`.  All arithmetic is integral and the verifier
rejects any bound above 5000; the certified cases stay below it.  Every
coefficient not assigned to `(P17F.6)` is checked nonnegative.  Pattern
seven changes only the label-four sign and passes the same construction
with its exact base coefficients.

**Proposition 17F (all eight-defect sign sectors).**  Under the
deficit-eight hypotheses, `g_(S-8)>=0` for every normalized sign pattern.

**Proof.**  Proposition 17C covers patterns zero through three,
Proposition 17D covers pattern four, Proposition 17E covers pattern five,
and `(P17F.1)--(P17F.8)` cover patterns six and seven.  The exact replay
reports a pass for every sector.  QED.

The consolidated replay is

```text
for p in 0 1 2 3; do
  ./analyze_su2_outer_eight 1 recursive $p
done
for p in 4 5 6 7; do
  ./analyze_su2_outer_eight 1 deficit8-chambers $p
done
```

and is recorded in `certificates/su2_outer_eight_all_patterns.log`.  Thus
the proved outer frontier is now every sign sector through deficit eight.
The next unresolved outer coefficient is deficit ten.

### Deficit-ten finite reduction

The same outer peeling identity now removes every factor of label at least
six: re-adjoining such a factor changes `g_(S-10)` by a sum of the proved
outer coefficients at deficits zero, two, four, six, and eight.  Hence the
deficit-ten problem reduces exactly to labels `1,2,3,4,5`.  Central-parity
normalization fixes the label-one sign and leaves sixteen sign patterns.

For counts `(c_1,...,c_5)`, admissibility is

```text
sum_(j=1)^5 j c_j - 10 >= max{j:c_j>0}.            (D10.1)
```

This set is upward closed and has 94 product-order minimal bases.  The
enumeration with total label at most twenty is exhaustive.  If a minimal
base had total label at least 21, deleting any present factor would leave
total label at least 16, still at least ten plus the largest possible
remaining label five, contradicting minimality.

The total count-degree of `g_(S-10)` is at most ten.  Indeed an invariant
selected subset of total label `2j<=10` contributes selection degree at
most `2j`, while the complementary near-top multiplicity has degree at most
`5-j`; their sum is `5+j<=10`.

The strict C++ analyzer
`character_ring_iter/analyze_su2_outer_ten.cpp` implements the exact subset
contraction with arbitrary-precision integers.  Its initial finite sweep
passed 846,112 admissible count/sign cases through twenty factors, with
minimum value `9` at counts `(0,0,0,0,3)` in the all-plus sector.  This is
diagnostic evidence only, not an unbounded deficit-ten certificate.  The
next task is to classify the Newton failures over the 94 minimal orthants
and find the analogue of the three ratio chambers used at deficit eight.

The direct Newton classification now proves four sectors uniformly.  For
each of patterns zero through three, every coefficient of the degree-ten
Newton polynomial is nonnegative on all 94 minimal orthants: there are
`94*C(15,5)=282282` exact coefficients per pattern.  These four patterns
are precisely those with positive normalized signs on labels two and
three; the signs on labels four and five are arbitrary.

**Proposition 17G (one quarter of the ten-defect sign sectors).**  Under
the deficit-ten hypotheses, `g_(S-10)>=0` whenever the normalized signs on
all label-two and label-three factors are positive.

**Proof.**  The total degree is at most ten, so the finite Newton expansion
on each minimal orthant is exact.  The strict verifier finds all 282,282
coefficients nonnegative in each of the four sign sectors.  The 94
orthants cover the admissible domain.  QED.

The remaining twelve sectors fail this direct certificate, all first at
the minimal base `(0,0,0,0,3)`.  Their first negative coefficients are

```text
patterns 4,5:    degree (1,1,1,1,0), value -1;
patterns 6,7:    degree (1,2,1,0,0), value -2;
patterns 8--11:  degree (1,1,2,0,0), value -1;
patterns 12--15: degree (4,1,0,0,1), value -1.     (D10.2)
```

These are failures of coefficientwise positivity, not negative values of
`g_(S-10)`.  The exact classification is recorded in
`certificates/su2_outer_ten_newton_classification.log`.  Recursive orthant
splitting and ratio-sensitive chambers are the next certificate layer.

The first failed sector has an exact recursive orthant certificate.  At a
node with base `b` and a negative Newton coefficient, choose a free count
coordinate `x` and partition

```text
{x>=b_x}={x=b_x} disjoint_union {x>=b_x+1}.        (D10.3)
```

The first child removes `x` from the free-coordinate mask; the second shifts
the base by one and retains the mask.  This is an exhaustive induction, and
a leaf with nonnegative Newton coefficients proves its whole orthant.

**Proposition 17H (a fifth ten-defect sign sector).**  Under the
deficit-ten hypotheses, `g_(S-10)>=0` for normalized sign pattern four,
`(s_1,s_2,s_3,s_4,s_5)=(+,+,-,+,+)`.

**Proof.**  The recursive exact verifier covers all 94 minimal orthants with
4,684 nodes, 2,389 nonnegative leaves, and 2,250 slice/tail splits.  Its
maximum depth is 50 and its largest visited base is `(11,15,44,4,3)`.
Every successful branch terminates at a nonnegative leaf; no depth or node
cap is accepted as a certificate.  Therefore repeated use of `(D10.3)`
proves the original admissible domain.  QED.

The replay is recorded in
`certificates/su2_outer_ten_recursive_pattern4.log`.  Pattern five differs
only in the normalized label-five sign and is the next recursive sector.

The pattern-five replay has the identical tree counts: 4,684 nodes, 2,389
nonnegative leaves, and 2,250 splits, again with maximum depth 50 and maximum
base `(11,15,44,4,3)`.

**Proposition 17I (the sixth ten-defect sign sector).**  Under the
deficit-ten hypotheses, `g_(S-10)>=0` for normalized sign pattern five,
`(s_1,s_2,s_3,s_4,s_5)=(+,+,-,+,-)`.

**Proof.**  Apply the exhaustive recursive partition `(D10.3)` using the
pattern-five exact coefficients.  Every branch terminates at a nonnegative
Newton leaf, with the counts stated above.  QED.

The replay is recorded in
`certificates/su2_outer_ten_recursive_pattern5.log`.  Thus six of the
sixteen normalized deficit-ten sectors are proved.

Pattern six also has an exact recursive certificate: 4,694 nodes, 2,394
nonnegative leaves, and 2,255 splits, with maximum depth 50 and maximum base
`(11,15,44,4,3)`.

**Proposition 17J (the seventh ten-defect sign sector).**  Under the
deficit-ten hypotheses, `g_(S-10)>=0` for normalized sign pattern six,
`(s_1,s_2,s_3,s_4,s_5)=(+,+,-,-,+)`.

**Proof.**  The exact recursive tree exhausts `(D10.3)` at every internal
node and terminates only at nonnegative Newton leaves.  QED.

The replay is recorded in
`certificates/su2_outer_ten_recursive_pattern6.log`.

Pattern seven has the same exact recursive tree counts as patterns four and
five: 4,684 nodes, 2,389 nonnegative leaves, and 2,250 splits, with maximum
depth 50 and maximum base `(11,15,44,4,3)`.

**Proposition 17K (one half of the ten-defect sign sectors).**  Under the
deficit-ten hypotheses, `g_(S-10)>=0` for normalized sign pattern seven,
`(s_1,s_2,s_3,s_4,s_5)=(+,+,-,-,-)`.  Consequently patterns zero through
seven, equivalently every normalized sector with positive label-two sign,
are proved.

**Proof.**  Apply the exhaustive recursive partition `(D10.3)`.  Every
branch terminates at a nonnegative Newton leaf; the verifier accepts neither
its depth cap nor its node cap as a certificate.  Together with Propositions
17G--17J this proves all eight sectors with positive label-two sign.  QED.

The replay is recorded in
`certificates/su2_outer_ten_recursive_pattern7.log`.

The remaining patterns eight through fifteen have negative label-two sign.
For pattern eight, the same recursive algorithm rejects at its 10,000-node
cap on the first base `(0,0,0,0,3)`, and a linear slope fan in `(c_1,c_2)`
also fails outside the first cone.  These are certificate failures, not
negative values.  They identify a different asymptotic scale: in the upper
cone `c_2=2c_1+v`, the bivariate Newton tail has competing terms along
`v` of order `c_1^2`, so a parabolic rather than linear ratio chamber is the
next exact target.

The upper-cone obstruction nevertheless has a strictly positive parabolic
leading face.  Write `c_1=u`, `c_2=2u+v` and hold the other counts fixed.
At the first minimal base the six Newton coefficients of parabolic weight
ten, for degrees `(10-2j,j,0,0,0)`, `0<=j<=5`, are

```text
42, -14, 6, -2, 2, 10.                            (D10.4)
```

After substituting `v=z u^2`, taking the leading monomial terms, and scaling
by `10!/42`, their polynomial is

```text
H(z)=1-30z+360z^2-1200z^3+3600z^4+7200z^5.        (D10.5)
```

This polynomial is strictly positive on the nonnegative real axis.  Indeed,

```text
H(z)=(1-15z)^2+15z^2 K(z),
K(z)=9-80z+240z^2+480z^3.
```

The only positive critical point of `K` is
`alpha=(sqrt(3)-1)/6`, where

```text
K(alpha)=(241-120sqrt(3))/9>0;
241^2-3*120^2=14881>0.                             (D10.6)
```

The strict `parabolic-tail` mode recomputes `(D10.4)` from the character
contraction and audits the coefficient decomposition and the final exact
radical comparison.  This proves coercivity of the two-coordinate
parabolic leading face.

The simultaneous leading face also has a closed positive form.  Give the
upper-cone variables weights `(1,2,2,2,2)`, write

```text
z=v/u^2,                 w=(c_3+c_4+c_5)/u^2,
```

and scale the weight-ten monomial face by `10!`.  Exact conversion of all
3,003 Newton coefficients gives

```text
F(z,w)=42 H(z)
 +1260 B(z)w+12600 C(z)w^2
 +50400(1+12z^2)w^3+75600(1+2z)w^4+30240w^5,     (D10.7)

B(z)=720z^4+120z^2-16z+1
    =(1-8z)^2+56z^2+720z^4,
C(z)=72z^3+36z^2-6z+1
    =(1-3z)^2+27z^2+72z^3.                       (D10.8)
```

Equations `(D10.5)--(D10.8)` show `F(z,w)>0` for every `z,w>=0`.
The full leading face depends on the three higher-label variables only
through their sum `w`; the `parabolic-leading-full` mode checks all 113
nonzero monomial coefficients against `(D10.7)`.  Thus simultaneous
parabolic growth of labels three through five is no longer an asymptotic
obstruction.  What remains is an exact domination of the lower-weight
monomials and the bounded part before that domination becomes effective.

That domination can be implemented without a large multidimensional box.
On each transformed minimal orthant, retain the first two chamber variables
`u,v` symbolically and expand the other three variables in the Newton basis:

```text
g_(S-10)=sum_d Q_d(u,v)
                    binom(x_3,d_3)binom(x_4,d_4)binom(x_5,d_5).  (D10.9)
```

There are 56 bivariate blocks per exceptional orthant.  For a block with a
negative monomial coefficient, give `v` weight one in the middle chamber
and weight two in the upper chamber.  If `W` is its maximum weighted degree,
write its leading face as `u^W L(v/u^s)`, with `s=1` or `2`.  Exact Sturm
sequences prove `L(z)>0` for `z>0`.  For every lower-weight negative monomial
`-C u^a v^b`, a second Sturm comparison supplies a dyadic bound

```text
L(z)>=2^r z^b.                                      (D10.10)
```

Consequently the leading face dominates all negative lower-weight terms
after an explicitly computed power-of-two threshold in `u`.  Below that
threshold, fixing the integer `u` leaves a univariate polynomial in `v`;
the verifier first accepts a nonnegative Newton expansion, and otherwise an
exact Sturm sequence proves it positive on the positive real axis (zeros at
the origin are removed as monomial factors).  Thus each `Q_d` is nonnegative
on the nonnegative integer quadrant, and `(D10.9)` proves its whole orthant.

For pattern eight, the lower chamber passes directly on all 78 transformed
minimal bases.  In the middle chamber, all 54 bases pass through 3,024
bivariate blocks, of which 421 require the Sturm/tail argument; the maximum
tail threshold is 256.  In the upper chamber, all 57 bases pass through
3,192 blocks, of which 570 require the argument; the maximum threshold is
1,024.

**Proposition 17L (a ninth ten-defect sign sector).**  Under the deficit-ten
hypotheses, `g_(S-10)>=0` for normalized sign pattern eight,
`(s_1,s_2,s_3,s_4,s_5)=(+,-,+,+,+)`.

**Proof.**  The lower, middle, and upper cones cover the complete
nonnegative `(c_1,c_2)` quadrant.  Their transformed minimal bases cover the
admissible count domain.  Direct Newton positivity proves the lower cone,
and `(D10.9)--(D10.10)` prove every middle and upper orthant by exact rational
Sturm arithmetic.  QED.

The replay is recorded in
`certificates/su2_outer_ten_bivariate_pattern8.log`.  Nine of the sixteen
normalized deficit-ten sectors are now proved.

Pattern nine has the identical chamber, block, and threshold counts.

**Proposition 17M (a tenth ten-defect sign sector).**  Under the deficit-ten
hypotheses, `g_(S-10)>=0` for normalized sign pattern nine,
`(s_1,s_2,s_3,s_4,s_5)=(+,-,+,+,-)`.

**Proof.**  Apply `(D10.9)--(D10.10)` with the pattern-nine coefficients.
Every bivariate block in every transformed minimal orthant passes.  QED.

The replay is recorded in
`certificates/su2_outer_ten_bivariate_pattern9.log`.  Ten of the sixteen
normalized deficit-ten sectors are now proved.

**Proposition 17N (an eleventh ten-defect sign sector).**  Under the
deficit-ten hypotheses, `g_(S-10)>=0` for normalized sign pattern ten,
`(s_1,s_2,s_3,s_4,s_5)=(+,-,+,-,+)`.

**Proof.**  The same three chamber cover and bivariate Sturm certificate
passes for every transformed minimal orthant.  QED.

The replay is recorded in
`certificates/su2_outer_ten_bivariate_pattern10.log`.

Pattern eleven has the same chamber and threshold counts as pattern ten.

**Proposition 17O (three quarters of the ten-defect sign sectors).**  Under
the deficit-ten hypotheses, `g_(S-10)>=0` for normalized sign pattern eleven,
`(s_1,s_2,s_3,s_4,s_5)=(+,-,+,-,-)`.  Consequently every normalized sector
with positive label-three sign is proved, irrespective of the signs on
labels two, four, and five.

**Proof.**  Apply the exact bivariate certificate on all three chambers.
Together with Propositions 17G--17N this covers patterns zero through eleven.
QED.

The replay is recorded in
`certificates/su2_outer_ten_bivariate_pattern11.log`.  Twelve of the sixteen
normalized deficit-ten sectors are now proved.

Pattern twelve requires the bivariate certificate on 13 of the 78 lower
minimal bases as well as on every middle and upper base.  All pass; the
maximum thresholds are respectively 64, 256, and 1,024.

**Proposition 17P (a thirteenth ten-defect sign sector).**  Under the
deficit-ten hypotheses, `g_(S-10)>=0` for normalized sign pattern twelve,
`(s_1,s_2,s_3,s_4,s_5)=(+,-,-,+,+)`.

**Proof.**  Apply `(D10.9)--(D10.10)` to each exceptional lower, middle, and
upper orthant.  Every exact Newton/Sturm block passes.  QED.

The replay is recorded in
`certificates/su2_outer_ten_bivariate_pattern12.log`.  Thirteen of the
sixteen normalized deficit-ten sectors are now proved.

**Proposition 17Q (a fourteenth ten-defect sign sector).**  Under the
deficit-ten hypotheses, `g_(S-10)>=0` for normalized sign pattern thirteen,
`(s_1,s_2,s_3,s_4,s_5)=(+,-,-,+,-)`.

**Proof.**  Pattern thirteen has the identical exact chamber, block, and
threshold counts as pattern twelve.  QED.

The replay is recorded in
`certificates/su2_outer_ten_bivariate_pattern13.log`.

**Proposition 17R (a fifteenth ten-defect sign sector).**  Under the
deficit-ten hypotheses, `g_(S-10)>=0` for normalized sign pattern fourteen,
`(s_1,s_2,s_3,s_4,s_5)=(+,-,-,-,+)`.

**Proof.**  The exact lower, middle, and upper bivariate certificates pass
on all transformed minimal bases.  QED.

The replay is recorded in
`certificates/su2_outer_ten_bivariate_pattern14.log`.

Pattern fifteen has the identical exact statistics as pattern fourteen.

**Proposition 17S (all ten-defect sign sectors).**  Under the deficit-ten
hypotheses, `g_(S-10)>=0` for every normalized sign pattern.

**Proof.**  Propositions 17G--17K prove patterns zero through seven,
Propositions 17L--17O prove patterns eight through eleven, and Propositions
17P--17R together with the pattern-fifteen replay prove patterns twelve
through fifteen.  For the last pattern the normalized signs are
`(s_1,s_2,s_3,s_4,s_5)=(+,-,-,-,-)`, and every lower, middle, and upper
bivariate block passes.  QED.

The pattern-fifteen replay is recorded in
`certificates/su2_outer_ten_bivariate_pattern15.log`.  The proved outer
frontier is now every sign sector through deficit ten.  The next unresolved
outer coefficient is deficit twelve.

All outer layers can be assembled into one palindromic polynomial.  Put

```text
P_q(t)=1+t+...+t^q
```

and retain the notation `q_i,s_i,S` of Proposition 16.  Define

```text
A_R(t)=sum_(T subset [n])
          (product_(i in T)s_i) N(T) t^(q(T)/2)
          product_(i notin T) P_(q_i)(t),               (P18.1)

q(T)=sum_(i in T)q_i.
```

Only invariant subsets occur, so every exponent `q(T)/2` in `(P18.1)` is an
integer.

**Proposition 18 (outer generating polynomial).**  If

```text
A_R(t)=sum_(d=0)^S a_d t^d,
```

then `A_R` is palindromic, `a_d=a_(S-d)`, and for every nonnegative target
`S-2d`,

```text
g_(S-2d)=[t^d](1-t)A_R(t)=a_d-a_(d-1),           (P18.2)
```

where `a_(-1)=0`.  Equivalently,

```text
A_R(t)=integral_SU(2) product_i[
          P_(q_i)(t)+s_i t^(q_i/2)chi_(q_i)(y)]dmu(y).  (P18.3)
```

Consequently `(OTM)` is exactly the initial-unimodality assertion

```text
a_0<=a_1<=...<=a_d
whenever S-2d>=max_i q_i.                         (OUM)
```

**Proof.**  Fix the subset `T` assigned to the integrated `y` coordinate.
Its contribution is its sign times `N(T)` times the product of the remaining
characters.  If the complementary labels have sum `L=S-q(T)`, their torus
weight polynomial, read downward from the highest weight, is

```text
product_(i notin T) P_(q_i)(t).
```

For any `SU(2)` representation, an irreducible multiplicity is the
difference between two adjacent weight multiplicities.  The coefficient of
`V_(S-2d)=V_(L-2[d-q(T)/2])` is therefore the coefficient of `t^d` in

```text
(1-t)t^(q(T)/2) product_(i notin T)P_(q_i)(t).
```

Summing `T` proves `(P18.2)`.  Expanding the product under the integral gives
`(P18.3)`.

Finally, each summand of `(P18.1)` has lowest degree `q(T)/2` and highest
degree

```text
q(T)/2+S-q(T)=S-q(T)/2;
```

the product of the `P_q` is symmetric between those endpoints.  Every
summand is therefore palindromic about degree `S/2`, proving the first
claim.  The equivalence with `(OUM)` follows from `(P18.2)` and the choice of
a largest removed factor in Proposition 15.  QED.

This formulation identifies the likely missing object: a graded complex or
Lefschetz module whose signed Hilbert polynomial is `A_R(t)` and whose
primitive Euler characteristics remain nonnegative through the outer range.
Ordinary coefficientwise positivity or monotonicity under adjoining a
factor is not sufficient, as the preceding exact failures show.

In fact the unrestricted partial contraction gives an equivalent global
form of the whole conjecture.

**Proposition 19 (partial-character and unimodality equivalence).**  The
central-character `Q3` theorem is equivalent to either of the following
statements, asserted for every signed remainder `R`:

```text
(PCP)  G_R(x)=integral R(x,y)dmu(y) is a genuine SU(2) character;

(FUM)  A_R(t) has nonnegative palindromic unimodal coefficients.
```

**Proof.**  Suppose `Q3` holds and write

```text
G_R=sum_(p>=0) g_p chi_p.
```

Let `epsilon` be the exchange parity of `R`.  For any `p>=1`, append the
factor `chi_p(x)+epsilon chi_p(y)`.  Its total number of minus factors is
even, and the exchange argument `(P15.4)` says that its `Q3` integral is
`2g_p`.  Hence every `g_p` is nonnegative.  If `epsilon=+1`, the same
argument with the plus trivial factor gives `g_0>=0`; if `epsilon=-1`, then
exchange antisymmetry gives `g_0=integral integral R=0`.  Thus `(PCP)`
follows.

Conversely, remove any factor from a nonempty even-minus signed list.  Its
sign equals the exchange parity of the remainder, so `(P15.4)` and `(PCP)`
make the original integral twice a nonnegative irreducible coefficient.
The empty case is immediate.  Hence `(PCP)` implies `Q3`.

Finally, Proposition 18 gives

```text
g_(S-2d)=a_d-a_(d-1)
```

from the outer coefficient down to the middle, while `A_R` is palindromic
and `a_0=1`.  Nonnegativity of all `g_p` is therefore equivalent to
nondecreasing coefficients up to the middle, which by palindromy is exactly
`(FUM)`.  Such a sequence is automatically nonnegative.  QED.

The strict exact C++ program
`character_ring_iter/search_su2_partial_character_cone.cpp` constructs
`G_R` directly in the signed two-coordinate representation ring and checks
every irreducible coefficient.  It first passed all 61,182 signed multisets
with labels at most five and at most seven factors, and the larger run passed
all 471,040 signed multisets with labels at most six and at most eight
factors.  This is evidence for
`(PCP)`, not its missing uniform proof.

There is an immediate unbounded sector in which `(PCP)` and `Q3` are
termwise positive.

**Proposition 20 (central-parity-aligned sector).**  Suppose every minus
factor has odd label and every plus factor has even label.  Then the signed
product satisfies `Q3`, with no bound on the labels or number of factors.
The same assertion holds uniformly in every finite fusion ring `SU(2)_k`.

**Proof.**  In the subset expansion of the signed integral, a subset `T`
assigned to one coordinate contributes

```text
(product_(i in T)s_i) N(T)N(T^c),
```

where `s_i=-1` precisely for a minus factor.  A necessary condition for
`N(T)` to be nonzero is that `sum_(i in T)q_i` be even.  Under the stated
parity alignment, the odd labels in `T` are exactly its minus factors.
Their number is therefore even, so `product_(i in T)s_i=+1`.  Every nonzero
summand is a product of two nonnegative invariant multiplicities.  QED.

At finite level the fusion rule has the same `Z/2` grading: every constituent
of `e_a e_b` has parity `a+b`.  Hence a subset with nonzero level-`k`
invariant multiplicity again has even total label, and the identical
termwise proof applies.

The same argument shows directly that the partial character in Proposition
19 is genuine for any remainder whose signs obey `s_i=(-1)^(q_i)`: every
invariant contracted subset has positive sign.

The partial-character equivalence is categorical and therefore has an exact
finite-level version.  Let `R_k` be the `SU(2)_k` fusion ring with basis
`e_0,...,e_k`, and let `tau_k` extract the coefficient of `e_0`.  For a
signed remainder define

```text
G_R^(k)=(id tensor tau_k)
          product_i[e_(q_i) tensor e_0
                    +s_i e_0 tensor e_(q_i)].           (P21.1)
```

**Proposition 21 (finite partial-character equivalence).**  The full
`GKS2*` condition in `SU(2)_k` is equivalent to

```text
G_R^(k) belongs to the nonnegative fusion cone
for every signed remainder R.                           (PCP_k)
```

**Proof.**  The fusion trace is Frobenius and the simple basis is
orthonormal:

```text
tau_k(e_p e_r)=indicator_(p=r).
```

If `epsilon` is the exchange parity of the remainder, append
`e_p tensor e_0+epsilon e_0 tensor e_p`.  Swapping the two tensor factors
shows exactly as in `(P15.4)` that the resulting signed trace is twice the
coefficient of `e_p` in `G_R^(k)`.  Thus `GKS2*` implies every coefficient
in `(PCP_k)` is nonnegative; the trivial coefficient is handled as in
Proposition 19.  Conversely, remove any factor from an even-minus signed
product.  Its sign equals the exchange parity of the remainder, so the same
identity and `(PCP_k)` give its signed trace as twice a nonnegative fusion
coefficient.  QED.

Consequently a uniform proof of `(PCP_k)` for all `k` proves the preferred
finite-ring theorem directly, without a packet-by-packet affine-wall split.
Corollary 3 then passes it exactly to ordinary `SU(2)` once the level exceeds
the explicit stabilization bound.

The strict exact C++ verifier
`character_ring_iter/search_su2_level_partial_character_cone.cpp` checked
every coefficient in `(PCP_k)` for all signed multisets with at most seven
factors at levels `1<=k<=6`.  It passed 235,514 cases in total, including
141,568 cases at level six.  This substantially strengthens the earlier
finite-level scalar checks, but it remains bounded evidence rather than a
uniform proof of `(PCP_k)`.

There is a uniform standard-monomial reduction of the kernel which a
positive complex must resolve.  Work first in the ordinary category and put
`V_i=Sym^(q_i)(C^2)`.  For a target `p` and a subset `T` of the signed
remainder positions, define

```text
I_T=Hom_(SL_2)(1, tensor_(i in T)V_i),
C_(p,T)=Hom_(SL_2)(V_p, tensor_(i notin T)V_i),
U_(p,T)=I_T tensor C_(p,T).                           (P22.1)
```

Tensoring the two intertwiners and restoring the original factor order gives
a canonical multiplication map

```text
mu_T: U_(p,T) -> C_p,
C_p=Hom_(SL_2)(V_p,tensor_i V_i).                    (P22.2)
```

For fixed `T`, this map identifies `U_(p,T)` with the channel in which the
`T` factors fuse through the tensor unit.  In particular `mu_empty` is the
identity of `C_p`.

Split the subsets according to the sign of their term in the partial
character:

```text
E_p=direct-sum_(|T intersect M| even) U_(p,T),
O_p=direct-sum_(|T intersect M| odd) U_(p,T),
mu_-=sum_(odd T) mu_T:O_p -> C_p.                    (P22.3)
```

**Proposition 22 (universal multiplication-kernel target).**  The target
coefficient of the partial character is

```text
g_p=dim E_p-dim O_p.                                 (P22.4)
```

Put

```text
kappa_p=dim ker(mu_-),
P_p=sum_(nonempty T, |T intersect M| even) dim U_(p,T).
```

Then

```text
kappa_p<=P_p                                         (UKS)
```

is sufficient for `g_p>=0`.  A proof of `(UKS)` for every signed list and
target proves `(PCP)` and hence `Q3`.

**Proof.**  The dimension of `U_(p,T)` is

```text
N(T) multiplicity_(V_p)(tensor_(i notin T)V_i),
```

so summing with its subset sign is precisely the coefficient formula for
`g_p`, proving `(P22.4)`.  Since the empty positive summand is `C_p`,

```text
dim E_p=dim C_p+P_p.
```

On the negative side,

```text
dim O_p=rank(mu_-)+kappa_p<=dim C_p+kappa_p.
```

Subtracting proves the assertion under `(UKS)`.  QED.

The preceding target is now also checked against the original Plucker
multiplication map, rather than only against a toric special fibre.  The
strict exact program
`character_ring_iter/analyze_su2_plucker_kernel.cpp` enumerates standard
two-row tableaux, straightens bracket graphs with

```text
p_(ad)p_(bc)=p_(ac)p_(bd)-p_(ab)p_(cd),  a<b<c<d,
```

and computes the rational rank of `mu_-`.  The support-disjoint sweep with
labels and targets at most four, at most seven remainder factors, and total
degree at most eighteen checked 3,201 cases.  Among the 640 source-deficit
cases, 474 had nonzero kernel; the maximum kernel was 231 and the minimum
slack `P_p-kappa_p` was zero.  Every case passed `(UKS)`.  This is bounded
evidence, not a proof.

It also changes the diagnosis of the first GT obstruction.  For

```text
[1^-,1^-,2^+,2^+,3^+],  target 1,
```

the original odd source has dimension eight and rank six, so its kernel has
dimension two, not the six-dimensional kernel of the fixed GT fibre.  The
positive channel on the two negative vertices has dimension two and pays
exactly these two relations.  Every pair of odd channel images has zero
intersection.  The relations are therefore genuinely higher-order.  At the
bracket level they are instances of

```text
p_(ac)p_(bd)p_(ef)-p_(ad)p_(bc)p_(ef)
-p_(ae)p_(bf)p_(cd)+p_(af)p_(be)p_(cd)=0.            (P22.5a)
```

Indeed `(P22.5a)` is the difference of the two Plucker quadrics on
`(ab,cd)` and `(ab,ef)`, after multiplication by the unused edge; their
common monomial `p_(ab)p_(cd)p_(ef)` cancels.  In the displayed obstruction
that common monomial lies in the two-dimensional positive channel.  Thus a
pair-intersection or XOR complex cannot be the final object: the next map
must contain these first Plucker syzygies.

The target-zero part of every odd-sign remainder is already exact, in both
the ordinary and finite fusion rings.

**Lemma 22H0 (odd-sign target-zero equality).**  If the signed remainder has
an odd number of minus factors, then at target `p=0`

```text
mu_-:O_0->C_0 is onto,              dim ker(mu_-)=P_0. (P22.5b0)
```

The same equality holds in every `SU(2)_k`.

**Proof.**  The subset containing every remainder factor is odd.  Its
channel is

```text
I_all tensor Hom(1,1)=C_0,
```

and its multiplication map is the identity, so `mu_-` is onto.  On the
other hand, the target-zero coefficient is the double fusion trace of the
signed remainder.  Exchanging the two tensor coordinates changes its sign
because the number of minus factors is odd, while leaving the trace fixed;
hence `g_0=0`.  Formula `(P22.4)` now gives

```text
0=dim C_0+P_0-dim O_0.
```

Surjectivity yields
`dim ker(mu_-)=dim O_0-dim C_0=P_0`.  The finite fusion trace is again
exchange-invariant and the all-factor channel is again the identity, so the
same proof applies at every level.  QED.

There is a rigorous relative presentation when one negative factor is
fundamental.  It isolates the exact higher-syzygy statement which remains.
Let `a` be the negative vertex, assume `d_a=1`, and let `R_d` be the
multidegree-`d` piece of the bracket algebra.  With the target vertex already
adjoined, set

```text
F_0(d)=direct-sum_(j!=a) R_(d-e_a-e_j) e_j,
F_1(d)=direct-sum_(i<j<l; a notin {i,j,l})
       R_(d-e_a-e_i-e_j-e_l) e_(ijl),                (P22.5b)
```

omitting summands with a negative coordinate, and define

```text
phi(h e_j)=p_(aj)h,
d_1(h e_(ijl))=p_(jl)h e_i-p_(il)h e_j+p_(ij)h e_l. (P22.5c)
```

The brackets in `(P22.5c)` are oriented.

**Lemma 22H (fundamental-strand relative exact sequence).**  The sequence

```text
F_1(d) --d_1--> F_0(d) --phi--> R_d
```

is exact at `F_0(d)`.  Removing the unique `a`-edge from every odd channel
defines a canonical map `iota_-:O_p->F_0(d)` with
`mu_-=phi iota_-`, and hence a short exact sequence

```text
0 -> ker(iota_-) -> ker(mu_-)
  -> image(iota_-) intersect image(d_1) -> 0.          (P22.5d)
```

**Proof.**  The identity `phi d_1=0` is the three-term Plucker relation on
`a,i,j,l`.  Conversely, the standard straightening proof for the ideal
generated by the row brackets `p_(aj)` reduces every relation among those
generators by a relation of `(P22.5c)`; equivalently, `(P22.5c)` is the
degreewise first presentation of that Plucker ideal.  This proves exactness
at `F_0(d)`.  Since `d_a=1`, every graph in `R_d` has a unique edge incident
with `a`; removing it also proves that `phi` is onto.  Every odd channel
contains `a`, and removing that unique edge from its graphs gives
`iota_-`, and restoring it gives `mu_-=phi iota_-`.  The map from
`ker(mu_-)` to the intersection in `(P22.5d)` is `iota_-`; its kernel is
`ker(iota_-)`, and it is onto by the definition of the intersection.  QED.

Thus the fundamental case splits exactly into an allocation kernel and a
relative Plucker carrier.  Positive channels exclude `a`; stripping the
unique `a`-edge from their complementary factor similarly gives
`iota_+:P_p->F_0(d)`.  Put

```text
A=image(iota_-),  B=image(iota_+),  K=ker(phi)=image(d_1),
alpha_p=dim ker(iota_-).                              (P22.5e)
```

Then `(P22.5d)` gives the exact identity

```text
kappa_p=alpha_p+dim(A intersect K).                   (P22.5f)
```

**Lemma 22H1 (two-sided carrier-capacity criterion).**  The stronger
inequality

```text
alpha_p+dim((A intersect K)+(B intersect K)) <= P_p   (P22.5g)
```

implies `(UKS)` in the fundamental strand.

**Proof.**  The inclusion `A intersect K` into the sum in `(P22.5g)`,
together with `(P22.5f)`, gives `kappa_p<=P_p`.  QED.

Criterion `(P22.5g)` records the positive internal carrier `B intersect K`
even though only `A intersect K` is needed in `(P22.5f)`.  This makes it a
strictly stronger but more symmetric target.  It has an exact residual and
section form.  Put

```text
K_-=A intersect K,                    K_+=B intersect K,
beta_p=dim ker(iota_+),               r_p=rank(mu_+),
c_p=dim(K_-/(K_- intersect K_+)).                     (P22.5h)
```

Rank-nullity applied first to `mu_+=phi iota_+` and then to the restriction
of `iota_+` to `ker(mu_+)` gives

```text
P_p=beta_p+dim K_++r_p.                              (P22.5i)
```

Since `dim(K_-+K_+)=dim K_++c_p`, criterion `(P22.5g)` is therefore exactly

```text
alpha_p+c_p<=beta_p+r_p.                             (P22.5j)
```

Thus the unmatched negative allocation and carrier demand is compared with
the positive allocation kernel and global positive image; the common
carrier `K_- intersect K_+` and the whole internal positive carrier cancel
before any injection is sought.

For a section `s:R_d->F_0(d)` of `phi`, define the defect

```text
D_s=(1-s phi)iota_+:P_p->K.                          (P22.5k)
```

**Lemma 22H2 (optimal-section criterion).**  There is a section `s` for
which

```text
K_- is contained in image(D_s),       dim ker(D_s)>=alpha_p              (P22.5l)
```

if and only if

```text
c_p<=r_p,                             alpha_p+c_p<=beta_p+r_p.           (P22.5m)
```

When these inequalities hold, `s` can be chosen so that

```text
image(D_s)=K_-+K_+,
rank(D_s)=dim(K_-+K_+),
dim ker(D_s)=P_p-dim(K_-+K_+).                        (P22.5n)
```

**Proof.**  Fix one section `s_0`.  Every other section is `s_0+t`, where
`t:R_d->K` is arbitrary, and

```text
D_(s_0+t)=D_(s_0)-t mu_+.
```

On `ker(mu_+)` every defect is therefore the fixed map `iota_+`; its image
is `K_+` and its kernel has dimension `beta_p`.  On a complement `W` of
`ker(mu_+)`, the map `mu_+` is an isomorphism onto its image, so varying `t`
allows the defect values on the `r_p` basis directions of `W` to be chosen
arbitrarily in `K`.  Covering `K_-` modulo the already fixed `K_+` is
possible exactly when `c_p<=r_p`.  Under this condition choose `c_p` values
which lift a basis of `K_-/(K_- intersect K_+)` and choose every remaining
value inside `K_-+K_+`.  This gives the image and rank in `(P22.5n)`.
Equation `(P22.5i)` then shows that its kernel has dimension at least
`alpha_p` exactly when the second inequality in `(P22.5m)` holds.

Conversely, every defect image contains `K_+`, and its quotient modulo
`K_+` is generated by at most `r_p` complement directions.  The first
inequality is therefore necessary for `(P22.5l)`.  Any image covering
`K_-` has rank at least `dim(K_-+K_+)`; rank-nullity and `(P22.5i)` make the
second inequality necessary.  QED.

The second inequality in `(P22.5m)` is exactly the two-sided capacity
criterion and hence already proves `(UKS)` by Lemma 22H1.  The first says
when that dimension proof can additionally be realized by one optimal
row-ideal section.

The analyzer now constructs this presentation exactly.  In the first
pair-disjoint carrier case

```text
[1^-,2^+,2^+,2^+,3^+],  target 2,
```

it finds `kappa_p=1`, `P_p=7`, lifts the nine-term odd relation to four
`F_1` terms, and uses second-syzygy freedom to replace it by five terms whose
chosen faces all lie in the positive stripped span.  The first obstruction
to the *face-only* normalization in the wider replay is

```text
[1^-,2^+,2^+,2^+,2^+,3^+],  target 2.                (P22.5o)
```

Here `dim C_p=30`, `dim O_p=44`, `rank(mu_-)=25`, `kappa_p=19`, and
`P_p=41`.  Pair intersections and the four-term relations `(P22.5a)` span
12 kernel directions.  Modulo the directly positive `F_1` faces the whole
kernel has quotient rank three, while those 12 relations account for only
two of the three directions.  Thus the provisional face differential covers
18 of the 19 directions.  Deterministic row-ideal sections repair this
defect collectively, and one perturbed section together with the pair and
four-term relations has exact rank 19.  More intrinsically, Lemma 22H2
constructs a single optimal section without a random search.  The analyzer
finds

```text
alpha_p=2,  dim(A intersect K)=17,
dim(B intersect K)=13,
dim((A intersect K)+(B intersect K))=28,
beta_p=1,  c_p=15,  r_p=27,
15<=27,  2+15<=1+27.                                  (P22.5p)
```

Thus both `(UKS)` and the stronger criterion `(P22.5g)` pass in this first
face-only obstruction.  The large slack `41-19=22` confirms that the missing
face was a defect of the provisional differential rather than an
obstruction to the kernel bound.

There is a useful negative diagnostic for a simpler construction.  For a
positive standard bracket graph, remove its unique fundamental edge.  When
the component of the other endpoint misses the target, adding that component
to the cut in two ways gives two odd factorizations of the same graph and
hence a canonical kernel relation.  In `(P22.5o)` this raw cut-boundary map
has rank six.  Moving the target and fundamental vertices through 49 tested
standard-monomial orders raises the best rank to twelve, and even the union
of all 49 images still has rank twelve.  Thus raw component switching sees
the lower Plucker relations but not the seven remaining higher directions;
the row-ideal section in Lemma 22H2 is genuinely stronger.

The next bounded replay isolates the difference between arbitrary section
perturbations and the optimal section.  At

```text
[1^-,2^+,2^+,2^+,3^+,4^+],  target 2,                (P22.5q)
```

one has `dim C_p=34`, `dim O_p=36`, `rank(mu_-)=27`, `kappa_p=9`, and
`P_p=25`.  All six previously used deterministic sections, even after the
lower pair relations are added, span only four kernel directions.  Here,
however,

```text
alpha_p=1,  dim K_-=8,  dim K_+=0,
beta_p=0,   c_p=8,      r_p=25,
8<=25,      1+8<=0+25.                               (P22.5r)
```

The constructive proof of Lemma 22H2 produces an exact defect of rank eight
and kernel dimension seventeen.  Its carrier intersections, together with
the one-dimensional allocation kernel, span all nine true kernel
directions, and every resulting relation vanishes under the original
Plucker multiplication matrix.  Thus the optimal construction removes a
genuine failure of the ad hoc section family rather than merely improving
its numerical rank.

The multidegrees identify where to look next.  If `U` is the span of the
vertices other than `a`, then the generators of `F_0` transform as `U` and
those of `F_1` as `exterior^3 U`.  The exact second-syzygy verifier
`verify_su2_plucker_second_syzygies.cpp` finds two tableau forms for the next
generators: one repeated anchor plus three other vertices, and five distinct
vertices.  These are precisely the two weight types in the Schur module
`S_(2,1,1,1) U`.  The face-only residual in `(P22.5o)` is a 12-term
alternating relation supported on seven active vertices, and the section
defect supplies it exactly.  A bounded exact replay through six factors,
labels and targets at most four, and total label degree at most 14 checks 71
fundamental-strand cases, including 23 source-deficit cases, and verifies
the recovered relation span, `(P22.5j)`, and the two inequalities in
`(P22.5m)` in every case.  It now also checks the deletion--contraction
kernel-dimension identity below in every source-deficit case.  The minimum
slacks in both inequalities in `(P22.5m)` are zero, so neither admits a
uniform positive-margin weakening.  The remaining uniform task is to prove the
residual inequality `(P22.5j)` in arbitrary multidegree; proving also
`c_p<=r_p` realizes it by the optimal differential of Lemma 22H2.  By Lemma
22H0 only positive targets remain in the odd-sign strand.  The exact commands,
bounded replay, and normalized 12-term face-only residual are recorded in
`certificates/su2_plucker_relative_kernel.log`.

**Target 22H3 (positive-target residual section inequalities).**  For every
fundamental negative vertex and every target `p>0`, prove

```text
c_p<=r_p,                      alpha_p+c_p<=beta_p+r_p. (P22.5s)
```

The first inequality constructs a defect covering the full negative carrier
modulo the shared carrier; the second supplies its allocation capacity.  By
Lemmas 22H1--22H2, `(P22.5s)` proves `(UKS)` in the complete fundamental
strand and realizes the proof by one positive-channel differential.

At the free bracket-graph level there is a canonical model for this
differential.  Fix a graph and a component with odd minus parity.  The odd
cuts are the even cuts translated by that component; sending every nonempty
even cut to the difference between its translated odd cut and the base odd
cut identifies the even-cut space with the kernel of the odd-cut summation
map.  The rank-twelve obstruction above shows that this map does not descend
naively through Plucker straightening.  Target 22H3 asks for its derived
descent through the row-ideal resolution: the higher Plucker homotopies must
land in `K_-+K_+`, while `(P22.5s)` says that the transferred map uses no
more than the original positive-channel dimension.  This is the precise
homological-perturbation statement left by the computations.

There is an exact target-recursion behind this derived contraction.  Let
`star` denote the target vertex and let `bar(d)` be obtained from `d` by
deleting the fundamental vertex `a` and replacing `d_star=p` by `p+1`.
Write

```text
bar(F)_0=direct-sum_(j!=star)
             R_(bar(d)-e_star-e_j) e_j,
bar(phi)(h e_j)=p_(star,j)h.                           (P22.5t)
```

**Lemma 22H3A (fundamental-target deletion--contraction).**  Substitution of
the target vector for the fundamental vector gives exact rows

```text
0 -> R_(d-e_a-e_star) --p_(a,star)--> R_d
  --m--> R_(bar(d)) -> 0,

0 -> R_(d-e_a-e_star)e_star -> F_0(d)
  --pi--> bar(F)_0 -> 0,                               (P22.5u)
```

intertwined by `phi` and `bar(phi)`.  The left vertical map is multiplication
by `p_(a,star)`, both right vertical maps are onto, and projection induces a
canonical isomorphism

```text
pi:ker(phi) -> ker(bar(phi)).                          (P22.5v)
```

**Proof.**  In the bracket algebra, the ideal of the diagonal on which the
two vector variables `a` and `star` coincide is the principal ideal generated
by their determinant `p_(a,star)`.  Taking the displayed multidegree gives
the first exact row.  The second row merely separates the `e_star` summand
of `F_0(d)`.  For `j!=a,star`, substitution sends

```text
p_(aj)h  to  p_(star,j)m(h),
```

which is `(P22.5t)`, while the `e_star` summand maps to zero.  The left
vertical map identifies that summand with `ker(m)`.  Finally `bar(phi)` is
onto because `bar(d)_star=p+1>0`: every bracket monomial in the target
multidegree has an edge incident with `star`, and removing one such edge
gives a preimage in `bar(F)_0`.  The snake lemma now gives `(P22.5v)`.  QED.

Representation-theoretically, the two outer terms in the first row are the
`V_(p-1)` and `V_(p+1)` branches of
`V_1 tensor V_p=V_(p-1) direct-sum V_(p+1)`.  Under `(P22.5v)`, negative
carriers become the crossing row faces `p_(star,j)` created by merging a
cut that separated `a` from `star`; positive carriers become the
noncrossing row faces from cuts that kept them together.  Thus Target 22H3
can be restated entirely in the syzygy space `ker(bar(phi))`, with no
fundamental vertex left.  This is the natural ordinary precursor of the
two-branch Temperley--Lieb recursion; at the affine wall the `p+1` branch is
simply absent.

The isomorphism `(P22.5v)` does not permit the lower branch in `(P22.5u)` to
be discarded from the source-capacity argument.  For a plus word `Q`, put

```text
A_(r,s)=sum_(S subset Q)
          mult_r(tensor_(i in S)V_(q_i))
          mult_s(tensor_(i notin S)V_(q_i)).             (P22.5w)
```

**Lemma 22H3B (exact two-branch fundamental transfer).**  The coefficient
of `V_p` in the partial character with one negative fundamental and plus
word `Q` is

```text
g_p=A_(0,p-1)+A_(0,p+1)-A_(1,p),                        (P22.5x)
```

with an invalid label interpreted as zero.  In `SU(2)_k` the exact formula
is

```text
g_p^(k)=sum_(s in p star_k 1) A_(0,s)^(k)-A_(1,p)^(k).  (P22.5xk)
```

Consequently the complete ordinary fundamental strand is equivalent to

```text
A_(1,p)<=A_(0,p-1)+A_(0,p+1),
```

and at the affine wall `p=k` the finite target becomes
`A_(1,k)^(k)<=A_(0,k-1)^(k)`.

**Proof.**  Expand the all-plus two-coordinate character as

```text
F_Q=sum_(r,s) A_(r,s) chi_r tensor chi_s.
```

The numbers in `(P22.5w)` are its coefficients, and complementing `S`
shows `A_(r,s)=A_(s,r)`.  Multiply by
`chi_1 tensor 1-1 tensor chi_1` and extract the coefficient of
`chi_p tensor chi_0`.  The first term contributes
`A_(p-1,0)+A_(p+1,0)`, while the second contributes `-A_(p,1)`.
Symmetry gives `(P22.5x)`.  The same argument in the finite fusion ring
replaces the two ordinary neighbors by `p star_k 1`, proving
`(P22.5xk)`.  QED.

The sharp affine-wall cases containing the simple current close exactly.

**Lemma 22H3C (odd-sign simple-current wall equality).**  In `SU(2)_k`, let
a signed remainder have odd exchange parity and suppose one of its plus
factors is the invertible label `k`.  Then its partial-character coefficient
at the wall is exactly zero.  In particular, for one negative fundamental,

```text
g_k^(k)=A_(0,k-1)^(k)-A_(1,k)^(k)=0.                  (P22.5z)
```

**Proof.**  Remove one plus occurrence of `k` and write the resulting signed
two-coordinate product as

```text
H=sum_(r,s) h_(r,s)e_r tensor e_s.
```

Odd exchange parity gives `h_(s,r)=-h_(r,s)`, so every `h_(r,r)` is zero.
Since `e_k e_r=e_(k-r)`, the coefficient of `e_k tensor e_0` after
reattaching `e_k tensor 1+1 tensor e_k` is exactly

```text
h_(0,0)+h_(k,k)=0.                                    (P22.5za)
```

For the final assertion, denote the coefficient array of the all-plus word
remaining after removal of `k` by `A`.  Adjoining the removed factor changes
the array to

```text
A'_(r,s)=A_(k-r,s)+A_(r,k-s).                          (P22.5zb)
```

Therefore

```text
A'_(1,k)=A_(k-1,k)+A_(1,0),
A'_(0,k-1)=A_(k,k-1)+A_(0,1).
```

Both pairs agree by the complement symmetry `A_(r,s)=A_(s,r)`.  Formula
`(P22.5xk)` now gives `(P22.5z)`.  QED.

This is uniform in the level and number of factors.  It removes every wall
case of odd exchange parity containing the affine simple current from the
remaining finite target.  In the fundamental strand, only plus words
supported in `2,...,k-1` can have positive wall slack.

There is an equally exact reduction when the simple current carries a minus
sign.

**Lemma 22H3D (negative simple-current wall reduction).**  Remove one
negative factor `e_k tensor 1-1 tensor e_k` from an odd-sign remainder and
write the resulting even-sign product as

```text
H=sum_(r,s) h_(r,s)e_r tensor e_s.
```

Then the wall coefficient of the original remainder is

```text
g_k^(k)=h_(0,0)-h_(k,k).                               (P22.5zc)
```

If `T` is the product of the corresponding signed two-coordinate fusion
operators and `P=N_k tensor N_k`, the same difference is

```text
h_(0,0)-h_(k,k)
 =1/2 <e_(0,0)-e_(k,k),T(e_(0,0)-e_(k,k))>.            (P22.5zd)
```

**Proof.**  Multiplication by the removed negative factor contributes
`h_(0,0)` from its first-coordinate term and `-h_(k,k)` from its
second-coordinate term, proving `(P22.5zc)`.  Every fusion matrix commutes
with `N_k`, so `T` commutes with `P`; all are self-adjoint and
`P e_(0,0)=e_(k,k)`.  Expanding the quadratic form in `(P22.5zd)` gives two
copies of `h_(0,0)` on the diagonal and two copies of `h_(k,k)` off the
diagonal.  QED.

Thus the remaining simple-current wall target is the single-vector
diagonal dominance

```text
h_(0,0)>=h_(k,k)                                      (SCW)
```

for every even-sign support-disjoint remainder which does not contain a
second label `k`.  This is weaker than positivity of `T` on the whole
`P=-1` eigenspace and retains exactly the vector needed by the wall
coefficient.  The strict C++ mode
`search_su2_level_partial_character_cone simple-current-wall` checked
266,627 such even-sign remainders through level eight and seven remaining
factors.  Among 93,365 cases with `h_(k,k)>0`, none was tight and the minimum
value of `h_(0,0)-h_(k,k)` was two.  The numerical value two is the smallest
possible strict gap, not evidence of a larger analytic margin: for every
nonempty even-sign remainder, complementing a factor assignment pairs two
equal contributions to each diagonal coefficient `h_(r,r)`, so every such
coefficient and their difference are even.  The substantive bounded
observation is the absence of an active equality.  This is bounded evidence
for `(SCW)`, not its uniform proof; its exact output is in
`certificates/su2_simple_current_wall.log`.

Stratifying the same run by the number of minus factors in `H` gives active
counts `2,310`, `33,790`, `47,992`, and `9,273` for respectively zero, two,
four, and six minus factors; every stratum has minimum strict gap two.  Thus
`(SCW)` is already nontrivial before signed cancellation.  Its first
uniform subtarget is the unsigned boundary-channel inequality

```text
sum_S m_0(S)m_0(S^c) >= sum_S m_k(S)m_k(S^c),          (SCW0)
```

where `m_r(S)` is the multiplicity of `e_r` in the product over `S`.
Higher even-minus sectors require a signed refinement of `(SCW0)`, but are
not the source of the first wall obstruction.

The unsigned target is an outer-corner specialization of the fusion-Hankel
defect from Proposition 6, and its update loses the negative Hankel term.

**Lemma 22H3E (simple-current boundary-column reduction).**  Let `C_Q` be
the coefficient matrix of a finite all-plus word, put

```text
d_Q=C_Q e_0,                 E_Q=H(d_Q)-C_Q,
```

and let the level be `k`.  Then

```text
E_Q(a,k)=C_Q(k-a,0)-C_Q(a,k),                           (P22.5ze)
E_Q(k,k)=C_Q(0,0)-C_Q(k,k).                             (P22.5zf)
```

If a further plus label `q` is adjoined, then

```text
E_(qQ)(a,k)
 =sum_(u in a star_k q)E_Q(u,k)
    +E_Q(a,k-q)-E_Q(k-a,q),                            (P22.5zg)

E_(qQ)(k,k)=2E_Q(k-q,k).                                (P22.5zh)
```

**Proof.**  The simple-current rule gives
`N_r(a,k)=indicator_(r=k-a)`.  Hence

```text
H(d_Q)_(a,k)=d_Q(k-a)=C_Q(k-a,0),
```

which proves `(P22.5ze)` and its corner specialization.  Proposition 6
gives

```text
E_(qQ)=N_qE_Q+E_QN_q-H(E_Qe_q).
```

At `(a,k)`, the first two terms are respectively the fusion-interval sum
and `E_Q(a,k-q)`.  The simple-current rule used above makes the Hankel term
`(E_Qe_q)_(k-a)=E_Q(k-a,q)`, proving the first recurrence.  At `(k,k)`,
each of the first two terms is `E_Q(k-q,k)`.  The last term is
`(E_Qe_q)_0=E_Q(0,q)=0`, because every fusion-Hankel defect has zero row and
column at zero.  This proves the corner specialization.  QED.

Consequently `(SCW0)` for every nonempty word follows from the strictly
sharper boundary-column target

```text
C_Q(a,k)<=C_Q(k-a,0),             0<=a<=k.             (SCB)
```

The induction is exact: choose any last factor `q` and use `(P22.5zh)` at
`a=k-q`.  Target `(SCB)` is the two-layer switching defect with one terminal
equal to the affine simple current.  It retains much less state than the
full finite defect matrix and supplies exactly the datum required by the
outer-corner update.  A proof still needs a recurrence or injection for the
whole displayed boundary column.  By `(P22.5zg)`, its exact one-step form is
the finite fusion-concavity inequality

```text
E_Q(k-a,q)
 <=sum_(u in a star_k q)E_Q(u,k)+E_Q(a,k-q).            (SCBI)
```

Proving `(SCBI)` along an ordered construction of every plus word preserves
`(SCB)` and hence proves `(SCW0)`.  The negative term is now one specified
interior defect entry; no unspecified wall correction remains.

The complete boundary column can be proved uniformly on the finite
minuscule ray.

**Lemma 22H3F (anti-diagonal reflection for fundamental plus words).**  In
`SU(2)_k`, if `Q` consists of any number of plus labels `1`, then `(SCB)`
holds for every `0<=a<=k`.

**Proof.**  The coefficient `C_Q(a,b)` counts walks in the square graph
`{0,...,k}^2`, starting at `(0,0)`, in which each fundamental plus factor
moves exactly one coordinate by one edge of the `A_(k+1)` fusion graph.
Every walk ending at `(a,k)` meets the anti-diagonal `x+y=k`.  At its first
meeting, reflect the remaining suffix by

```text
R(x,y)=(k-y,k-x).
```

The map `R` fixes the anti-diagonal pointwise and is an automorphism of the
square graph, so the reflected suffix is again a valid walk with the same
number and order of steps.  Its endpoint is `R(a,k)=(0,k-a)`.  Reflection
after the first meeting is injective (indeed reversible on the image), and
coordinate exchange gives

```text
C_Q(a,k)<=C_Q(0,k-a)=C_Q(k-a,0).
```

This is `(SCB)`.  QED.

The preceding reflection is the boundary shadow of an exact finite
type-`C_2` transform.  It gives a stronger result and, crucially, survives
one arbitrary irreducible plus block.

**Lemma 22H3G (finite `C_2` transform and affine top-row cancellation).**
Fix a level `k`, put

```text
theta_j=(j+1)pi/(k+2),       x_j=2cos(theta_j),
w_j=2 sin(theta_j)^2/(k+2),  0<=j<=k,
mu_k=sum_j w_j delta_(x_j),
dnu_k(x,y)=(1/2)(x-y)^2 dmu_k(x)dmu_k(y).
```

For `0<=s<=r<=k-1`, let

```text
Psi_(r,s)(x,y)
 =[chi_(r+1)(x)chi_s(y)-chi_s(x)chi_(r+1)(y)]/(x-y),
D_a=Psi_(a-1,0).
```

Then the `Psi_(r,s)` form an orthonormal basis on the unordered
off-diagonal spectral grid.  For every finite plus word `Q` and
`1<=a,b<=k`,

```text
E_Q(a,b)=<D_a D_b F_Q>_(nu_k),
F_Q=product_(q in Q)[chi_q(x)+chi_q(y)].              (P22.5zj)
```

Writing `ell=k-1`, multiplication by the fundamental plus generator is the
nonnegative triangular-alcove adjacency

```text
Psi_(1,0)Psi_(r,s)
 =sum Psi_(u,v),                                      (P22.5zk)
```

where `(u,v)` runs through `(r+1,s),(r-1,s),(r,s+1),
(r,s-1)` which satisfy `0<=v<=u<=ell`.  At the affine top row one has the
two-term cancellation

```text
D_k[chi_q(x)+chi_q(y)]
 =Psi_(ell,q)+Psi_(ell-q,0),       1<=q<k.             (P22.5zl)
```

Consequently:

```text
E_[1^n](a,b)>=0                         (1<=a,b<=k),
E_[1^n,q](a,k)>=0                       (1<=a<=k,
                                         n>=0, 1<=q<k). (P22.5zm)
```

Thus the complete finite two-minus defect is nonnegative on every
fundamental plus ray, and the full simple-current boundary column remains
nonnegative after adjoining one arbitrary irreducible plus label.

**Proof.**  The finite fusion ring is the function algebra on the displayed
`SU(2)_k` spectral nodes, with the `chi_a` orthonormal for `mu_k`.  If
`C_Q` is the plus coefficient matrix and `d_Q=C_Qe_0`, then

```text
H(d_Q)_(a,b)
 =integral chi_a(x)chi_b(x)F_Q(x,y)dmu_k(x)dmu_k(y).
```

Subtract the analogous formula for `C_Q(a,b)` and average with its image
under `x<->y`.  This gives

```text
E_Q(a,b)
 =(1/2)integral [chi_a(x)-chi_a(y)]
                  [chi_b(x)-chi_b(y)]F_Q dmu_k(x)dmu_k(y),
```

which is `(P22.5zj)`.

After multiplication by `(x-y)`, the numerator defining `Psi_(r,s)` is the
exterior product of the two orthonormal vectors `chi_(r+1)` and `chi_s`.
Exterior-product orthogonality therefore makes the displayed `Psi_(r,s)`
orthonormal for `nu_k`.  There are

```text
k(k+1)/2
```

of them, exactly the number of unordered pairs of distinct spectral nodes,
so they form a basis.  The ordinary minuscule `USp(4)` rule gives the four
neighbors in `(P22.5zk)`.  At `r=ell`, its outward term is zero on the
finite grid because `chi_(k+1)(x_j)=0`; the remaining lower-wall terms are
deleted by dominance.  Hence `(P22.5zk)` is an exact nonnegative adjacency
rule on the finite triangle.

Finally `chi_k chi_q=chi_(k-q)` at level `k`, and hence

```text
(x-y)D_k[chi_q(x)+chi_q(y)]
 =chi_(k-q)(x)-chi_(k-q)(y)
   +chi_k(x)chi_q(y)-chi_q(x)chi_k(y).
```

The two terms after division by `x-y` are respectively
`Psi_(ell-q,0)` and `Psi_(ell,q)`, proving `(P22.5zl)`.  For a fundamental
word, `(P22.5zj)` is the coefficient of `D_a` in the nonnegative walk
`D_b Psi_(1,0)^n`, proving the first line of `(P22.5zm)`.  After adjoining
`q`, move `D_k[chi_q(x)+chi_q(y)]` together and use `(P22.5zl)`; the desired
coefficient is the sum of two nonnegative triangular-walk counts.  This
proves the second line.  QED.

The strict C++ mode `finite-c2-fundamental` compares `(P22.5zj)` against the
triangular adjacency through level thirty-two and one hundred fundamental
factors (1,155,440 exact entries).  The mode `finite-c2-wall-generator`
checks that `(P22.5zl)` has exactly its two stated unit coefficients through
level forty.  Their outputs are in
`certificates/su2_finite_c2_boundary.log`.

The tempting extension from `(P22.5zl)` to an invariant character-positive
wall cone is false.  Already at level three,

```text
D_3[chi_2(x)+chi_2(y)]^2
```

has coefficient `-2` at `Psi_(1,1)`, although all its row coefficients
needed for `(SCB)` are nonnegative.  The exact mode `finite-c2-wall-cone`
records this first obstruction.  A general proof must therefore preserve a
larger signed packet while controlling its row projection; ordinary finite
`C_2` character positivity cannot by itself close the induction.

The entire signed packet nevertheless has a closed allocation formula.  It
pinpoints the exact global switch which remains.

**Lemma 22H3H (affine cut-determinant expansion).**  For a subset `S` of a
finite plus word `Q`, let `m_S(r)` be the level-`k` fusion multiplicity of
`e_r` in the product over `S`, and put `T=Q\S`.  With `ell=k-1`,

```text
D_k F_Q
 =sum_(S,r,s: k-r!=s) m_S(r)m_T(s) sign(k-r-s)
    Psi_(max(k-r,s)-1,min(k-r,s)).                    (P22.5zo)
```

In particular its row coefficient is

```text
[Psi_(a-1,0)]D_kF_Q
 =sum_S [m_S(k-a)m_T(0)-m_S(k)m_T(a)]
 =E_Q(a,k).                                          (P22.5zp)
```

Thus negative `C_2` packets are exactly the allocation channels beyond the
affine cut `r+s=k`.  Boundary domination asks for a global switch from the
negative cut data `(r,s)=(k,a)` to the positive cut data
`(r,s)=(k-a,0)`.

**Proof.**  Expand the symmetric plus product as

```text
F_Q=sum_S A_S(x)A_T(y),
A_S=sum_r m_S(r)chi_r.
```

In the term containing `-chi_k(y)`, replace `S` by its complement.  Since
the simple current satisfies `chi_k chi_r=chi_(k-r)`, this gives

```text
D_kF_Q
 =sum_(S,r,s) m_S(r)m_T(s)
   [chi_(k-r)(x)chi_s(y)-chi_s(x)chi_(k-r)(y)]/(x-y).
```

If `k-r>s`, the quotient is `Psi_(k-r-1,s)`; if `k-r<s`, it is
`-Psi_(s-1,k-r)`; and it vanishes at equality.  This proves `(P22.5zo)`.
A row index has second coordinate zero.  Its positive occurrences have
`s=0,k-r=a`, while its negative occurrences have `r=k,s=a`, which proves
`(P22.5zp)`.  QED.

For one fixed cut, fusing a negative pair of trees `S->k`, `T->a` gives a
tree on all factors with output `k-a`, but different cuts need not remain
linearly independent after this reassociation.  The allocation and shell
counterexamples below show that this loss is real.  Hence `(P22.5zp)` calls
for a differential across the Boolean lattice of cuts, not independent
associativity injections for each cut.  This is the finite fusion analogue
of the Plucker row-ideal contraction in Target 22H3.

There is an equivalent exterior formulation which discards every packet
direction not seen by the desired row.

**Lemma 22H3I (simple-current exterior minor).**  Put

```text
K_Q(z)=product_(i in Q)(I+z_i N_(q_i))
```

and let `dGamma(N)=N tensor I+I tensor N` act on the exterior square of the
level-`k` fusion module.  Then

```text
E_Q(a,k)
 =[e_a wedge e_0]
   product_(i in Q)dGamma(N_(q_i))(e_k wedge e_0)       (P22.5zs)

 =[product_i z_i]
   det K_Q(z)[rows (a,0), columns (k,0)].               (P22.5zt)
```

The first update from the simple-current source is the positive two-packet

```text
dGamma(N_q)(e_k wedge e_0)
 =e_(k-q) wedge e_0+e_k wedge e_q.                     (P22.5zu)
```

If `J=N_k`, the source, target, and their entire orbit lie in the
`J wedge J=-1` sector.  Thus `(SCB)` is exactly positivity of the boundary
matrix coefficients in the cyclic exterior module generated by
`e_k wedge e_0`; no inequality on the complementary `C_2` packet directions
is required.

**Proof.**  Exterior powers turn `K_Q(z)` into the ordered product of the
exterior powers of its factors.  The coefficient linear in `z_i` in
`wedge^2(I+z_iN_(q_i))` is `dGamma(N_(q_i))`; hence selecting every variable
once proves the equality of the two displays.  The indicated minor is

```text
K_Q(a,k)K_Q(0,0)-K_Q(a,0)K_Q(0,k).
```

Since every fusion matrix commutes with the simple current,
`K_Q(a,k)=K_Q(k-a,0)`.  Its square-free coefficient is therefore

```text
sum_S [m_S(k-a)m_(S^c)(0)-m_S(a)m_(S^c)(k)],
```

which is `E_Q(a,k)` by `(P22.5ze)` and equals `(P22.5zp)` after complementing
the subsets in the negative term.  Finally
`N_qe_k=e_(k-q)` and `N_qe_0=e_q`, giving `(P22.5zu)`.  Commutation with
`J wedge J` proves the sector statement.  QED.

Formula `(P22.5zs)` also explains both the success and the limitation of
the affine top-row cancellation.  The first factor creates two positive
exterior carriers, but a later long fusion jump can cross the stationary
carrier and acquire the wedge sign.  The required Boolean-cut differential
is precisely a cancellation of those crossing histories inside this cyclic
sector.  The global sign-gauge obstruction below concerns the full exterior
square and therefore does not by itself decide this smaller boundary orbit.

The strict C++ mode `simple-current-exterior` constructs the additive
compounds directly and checks `(P22.5zs)` independently against the
fusion-Hankel defect.  Through level ten and seven plus factors it verifies
218,789 boundary coefficients and the `J wedge J=-1` symmetry in 1,113,397
exterior entries; every boundary coefficient is nonnegative.  The cyclic
sector is not entrywise positive: its first negative interior value is

```text
k=3, Q=[2,2], state=e_2 wedge e_1: value=-2.           (P22.5zv)
```

There are 165,463 negative interior entries in the full bounded run.  Thus
the evidence isolates a genuine boundary-positivity phenomenon rather than
an accidentally nonnegative transfer matrix.  The exact output is in
`certificates/su2_finite_c2_boundary.log`.

At odd level the cyclic exterior sector has a further exact factorization.
It shows that the wall problem is itself a full signed-correlation problem
for the simple-current orbit ring.

**Lemma 22H3J (odd-level orbit-ring equivalence).**  Let `k=2m-1`, let
`J=N_k`, and decompose the real fusion module as `V=V_+ direct-sum V_-`
under `J`.  For `0<=r<m`, use the orthonormal vectors

```text
u_r=(e_r+e_(k-r))/sqrt(2),
w_r=(-1)^r(e_r-e_(k-r))/sqrt(2).
```

The parity operator identifies `u_r` with `w_r`.  If `A_q` is the matrix of
`N_q` on `V_+` in the `u` basis, then `A_q` is entrywise nonnegative,
`A_(k-q)=A_q`, and on `V_-` one has

```text
N_q=(-1)^q A_q.                                      (P22.5zw)
```

Under the isometry

```text
Phi:V_+ tensor V_+ -> (wedge^2 V)_(J wedge J=-1),
Phi(u_r tensor u_s)=u_r wedge w_s,
```

the exterior transfer becomes

```text
Phi^(-1)dGamma(N_q)Phi
 =A_q tensor I+(-1)^q I tensor A_q,                  (P22.5zx)
```

and `Phi(u_0 tensor u_0)=e_k wedge e_0`.

For a word `Q`, put `epsilon_Q=sum_(q in Q)q mod 2`.  For each `0<=r<m`,
exactly one of the paired boundary targets `r,k-r` can be nonzero.  If

```text
a_r=r       when epsilon_Q == k-r mod 2,
a_r=k-r     otherwise,
```

then

```text
E_Q(a_r,k)
 =[u_r tensor u_0]
   product_(q in Q)[A_q tensor I+(-1)^q I tensor A_q]
   (u_0 tensor u_0),
E_Q(k-a_r,k)=0.                                      (P22.5zy)
```

The distinct matrices `A_0,...,A_(m-1)` form the nonnegative based orbit
fusion ring `O_k`; folding `q` and `k-q` gives its product coefficients.
Since `A_(k-q)=A_q` while `q` and `k-q` have opposite parity, either sign of
every nontrivial orbit generator can be realized in `(P22.5zx)`.  Therefore
the corner inequality `(SCW0)` for all plus words in `SU(2)_k` is
equivalent, when `k` is odd, to the full `GKS2*` property for `O_k`.  The
whole column `(SCB)` is equivalent to the strictly stronger assertion that
every partial-character coefficient of every signed orbit-ring word is
nonnegative.

**Proof.**  Let `P e_r=(-1)^r e_r`.  Fusion parity gives
`P N_q P=(-1)^qN_q`, while oddness of `k` gives `PJ=-JP`.  Thus `P` maps
`V_+` isometrically onto `V_-` and proves `(P22.5zw)`.  Acting on an orbit
sum `u_r`, the coefficient of `u_s` is the sum of the two nonnegative fusion
coefficients to `s` and `k-s`; hence every `A_q` is nonnegative.  Commutation
with `J` and the simple-current identity `N_(k-q)=JN_q` give
`A_(k-q)=A_q`.

The `J wedge J=-1` space is exactly `V_+ wedge V_-`, so `Phi` is an
isometry.  Applying `N_q` to its two wedge factors and using `(P22.5zw)`
gives `(P22.5zx)`.  The displayed formulas for `u_0,w_0` give
`e_k wedge e_0=u_0 wedge w_0`.

The ordered product in `(P22.5zx)` has coordinate-swap parity
`(-1)^(epsilon_Q)`.  Expanding `e_r wedge e_0` and
`e_(k-r) wedge e_0` into their `V_+ wedge V_-` parts shows that one target
selects the coefficient of `u_r tensor u_0`, while the other has the
opposite swap parity and vanishes.  This proves `(P22.5zy)`.  Restriction of
the fusion product to `V_+` proves closure and nonnegative structure
constants for the orbit ring.  Finally the two lifts `q` and `k-q` realize
the two signs.  At `r=0`, `(P22.5zy)` is the signed orbit-ring integral and
the parity-opposite target is zero, proving the `(SCW0)`--`GKS2*`
equivalence.  Allowing every `r` gives precisely partial-character
positivity and proves the final column statement.  QED.

The `simple-current-exterior` verifier independently checks `(P22.5zy)` in
40,101 odd-level word/target entries through level nine and seven factors.
This factorization is a structural reduction, not yet a proof: the orbit
rings are smaller, but they are non-group fusion rings (level three already
produces the Fibonacci ring), so their full `GKS2*` property still requires
an argument.

The orbit system in Lemma 22H3J is an explicit one-dimensional polynomial
hypergroup.

**Lemma 22H3K (tadpole spectral normal form).**  For `k=2m-1`, the
fundamental orbit matrix `A=A_1` acts on `u_0,...,u_(m-1)` by

```text
A u_0=u_1,
A u_r=u_(r-1)+u_(r+1)                 (1<=r<=m-2),
A u_(m-1)=u_(m-2)+u_(m-1).                         (P22.5zz)
```

Thus it is the path adjacency with one loop at the folded endpoint.  Put

```text
theta_j=(2j+1)pi/(2m+1),       x_j=2cos(theta_j),
omega_j=4 sin(theta_j)^2/(2m+1),       0<=j<m.       (P22.5zz1)
```

The `x_j` are the complete spectrum of `A`, the `omega_j` sum to one, and
the orbit characters are

```text
varphi_r(x_j)=chi_r(x_j),              0<=r<m.
```

They form an orthonormal real `GKS1` system for the measure
`sum_j omega_j delta_(x_j)`.  Consequently the full odd-level
simple-current boundary column is equivalent to the explicit inequalities

```text
sum_(i,j) omega_i omega_j varphi_r(x_i)
  product_t [varphi_(p_t)(x_i)
             +sigma_t varphi_(p_t)(x_j)] >=0,         (P22.5zz2)
```

for every `r`, every word `p_t in {1,...,m-1}`, and arbitrary signs
`sigma_t in {+1,-1}`.  The `r=0` specialization is exactly the corner
inequality `(SCW0)`, equivalently full `GKS2*` for the orbit ring; allowing
all `r` is the stronger partial-character positivity required by the whole
column `(SCB)`.

**Proof.**  Fold the ordinary `A_(2m)` fundamental fusion graph by its
endpoint reversal.  Every nonterminal orbit has the two path neighbors,
whereas the two middle vertices form one orbit and produce the terminal
loop, proving `(P22.5zz)`.  An eigenvector with first coordinate one has
coordinates `chi_r(2cos theta)`.  The loop boundary condition is

```text
chi_m(2cos theta)=chi_(m-1)(2cos theta),
```

whose roots are exactly the angles in `(P22.5zz1)`.  Folding the
`SU(2)_(2m-1)` sine-transform weight gives `omega_j`; alternatively the
spectral theorem at the cyclic vector `u_0` gives the same weights and their
sum one.  Since `chi_r(A)u_0=u_r`, spectral orthogonality proves that the
`varphi_r` are orthonormal.  Their nonnegative linearization is the folded
fusion product.  Finally `(P22.5zz2)` is `(P22.5zy)` written in this spectral
basis, with the two lifts `p_t,2m-1-p_t` realizing the sign `sigma_t`.
QED.

There is a second spectral model in which the folded fusion rule becomes
literal interval convolution on an odd cyclic group.  Unlike the tadpole
coordinates, it also turns the desired coefficient into one local discrete
mixed derivative.

**Lemma 22H3L (odd cyclic-gradient and Monge model).**  Keep `k=2m-1` and
put `n=k+2=2m+1`.  Relabel the orbit basis by its unique even lift:

```text
B_a=A_(min(2a,k-2a)),                    0<=a<m.
```

Then

```text
B_a B_b=sum_(c=|a-b|)^(min(a+b,2m-1-a-b)) B_c.       (P22.5zz3)
```

Let `C_n` be the group of `n`th roots of unity, write

```text
beta_a(u)=sum_(t=-a)^a u^t,
h(u)=2-u-u^(-1)=(1-u)(1-u^(-1)),
```

and identify `u` with `u^(-1)`.  On the `m` nontrivial inverse pairs put

```text
rho_m([u])=h(u)/n.                                  (P22.5zz4)
```

This is a probability measure, and `B_a` is represented by `beta_a`.
Equivalently, for every inversion-invariant function,

```text
<f>_(rho_m)=(1/(2n)) sum_(u^n=1) h(u)f(u).           (P22.5zz5)
```

The functions `beta_0,...,beta_(m-1)` are orthonormal for this measure.
Their cyclic gradients are the two-point roots

```text
(1-u)beta_a(u)=u^(-a)-u^(a+1).                       (P22.5zz6)
```

Now take any signed orbit word and set

```text
P(u,v)=product_i[beta_(a_i)(u)+sigma_i beta_(a_i)(v)]
      =sum_(r,s in Z/nZ) P_hat(r,s)u^r v^s.          (P22.5zz7)
```

If `g_c` is its partial-character coefficient at `B_c`, then

```text
g_c=P_hat(c,0)-P_hat(c,1)
       -P_hat(c+1,0)+P_hat(c+1,1),       0<=c<m.     (P22.5zz8)
```

Thus the corner `GKS2*` inequality is exactly

```text
P_hat(0,0)+P_hat(1,1)>=2P_hat(1,0),                  (P22.5zz9)
```

and the whole column `(SCB)` is the adjacent Monge family `(P22.5zz8)`.
The coefficient array in `(P22.5zz7)` is obtained from the unit mass at
`(0,0)` by the signed interval-convolution updates

```text
C -> (1_[-a,a] *_1 C)+sigma(1_[-a,a] *_2 C)          (P22.5zz10)
```

on `C_n times C_n`.

Finally, because `n` is odd, the change of variables

```text
u=zw,                     v=z/w
```

is a bijection of `C_n times C_n`.  In these coordinates every factor and
the Weyl weight have the rank-one expansions

```text
beta_a(zw)+beta_a(z/w)
 =2+sum_(t=1)^a(z^t+z^(-t))(w^t+w^(-t)),

beta_a(zw)-beta_a(z/w)
 =sum_(t=1)^a(z^t-z^(-t))(w^t-w^(-t)),

h(zw)h(z/w)
 =[(z+z^(-1))-(w+w^(-1))]^2.                        (P22.5zz11)
```

**Proof.**  The orbit of every label contains exactly one even label.
Restrict the level-`k` fusion rule for `2a` and `2b` to its even outputs
`2c`.  Dividing the ordinary lower and affine upper bounds by two gives
`(P22.5zz3)`.

The nontrivial elements of `C_n` are `m` inverse pairs.  Since

```text
sum_(u^n=1) h(u)=2n
```

and `h(1)=0`, `(P22.5zz4)` has total mass one and `(P22.5zz5)` follows.
Identity `(P22.5zz6)` is the telescoping sum of the interval polynomial.
For `0<=a,b<m`, it gives

```text
(1-u)beta_a(u)(1-u^(-1))beta_b(u)
 =(u^(-a)-u^(a+1))(u^b-u^(-b-1)).
```

The four exponents have absolute value below `n`; their cyclic constant
term is two when `a=b` and zero otherwise.  The factor `1/2` in
`(P22.5zz5)` therefore proves orthonormality.

Put `t_c(u)=u^c+u^(-c)`, including `t_0=2`.  A second telescoping identity
is

```text
h(u)beta_c(u)=t_c(u)-t_(c+1)(u).                    (P22.5zz12)
```

The polynomial `P` is invariant under inversion in either variable, so
uniform cyclic averaging gives

```text
(1/n^2)sum_(u,v) t_r(u)t_s(v)P(u,v)=4P_hat(r,s).
```

Apply `(P22.5zz5)` in both variables, insert `(P22.5zz12)` for
`h(u)beta_c(u)` and `h(v)=t_0(v)-t_1(v)`, and divide the last display by
four.  The result is exactly `(P22.5zz8)`.  Its `c=0` case is
`(P22.5zz9)`, while multiplication by an interval polynomial in either
variable is `(P22.5zz10)`.

The determinant of the exponent change `(z,w)->(zw,z/w)` is `-2`, which
is invertible modulo odd `n`; hence it is bijective.  Pairing the `t` and
`-t` terms in each interval polynomial proves the first two identities in
`(P22.5zz11)`.  The last is the elementary identity

```text
(2-zw-z^(-1)w^(-1))(2-z/w-z^(-1)w)
 =[(z+z^(-1))-(w+w^(-1))]^2.
```

This completes the proof.  QED.

Lemma 22H3L replaces the non-group tadpole multiplication by a group
convolution problem with only one residual local inequality.  It does not
by itself prove that inequality: after the change in `(P22.5zz11)`, plus
factors are positive-semidefinite rank-one kernels and minus factors are
negative-semidefinite sine kernels, so an even-minus product is
positive-semidefinite by the Schur product theorem; however the squared
Weyl difference is not a positive-semidefinite kernel (its diagonal is zero
but it is not the zero kernel).  The missing step is
therefore a precise adjacent-Monge theorem for products of the nested cyclic
interval kernels, rather than a general positive-semidefinite argument.

The strict C++ program `search_su2_odd_orbit_gks.cpp` independently evolves
both sides of `(P22.5zz8)`.  Its `cyclic-check` mode compares every orbit
boundary coefficient with the cyclic Monge difference.  The larger direct
orbit run exhausts support-disjoint sign choices; its exact output is
recorded in `certificates/su2_odd_orbit_gks.log`: ranks two through eight,
up to twelve factors, give 1,866,279 signed words and 14,328,179 boundary
coefficients, all nonnegative.  The independent cyclic model agrees in all
100,484 boundary entries through rank six and eight factors.  These are
bounded checks, not the missing adjacent-Monge theorem.

The rank-one expansion also has a finite three-charge form parallel to
Propositions 13--15.

**Lemma 22H3M (cyclic three-charge Gram and outer Turan reduction).**  For a
signed word in the `B_a` basis, make one local feature choice per factor:

```text
plus B_a:   either 1 with weight 2, or C_t(z)=z^t+z^(-t), 1<=t<=a;
minus B_a:  S_t(z)=z^t-z^(-t),                         1<=t<=a.
```

A complete choice `gamma` has weight `w_gamma` and cyclic Laurent
polynomial

```text
A_gamma(z)=sum_(r in Z/nZ) a_r(gamma)z^r.
```

If the number of minus factors is even, then the corner coefficient is

```text
g_0=sum_gamma w_gamma[
       a_0(gamma)^2-2a_1(gamma)^2
                         +a_0(gamma)a_2(gamma)].       (P22.5zz13)
```

Equivalently, there are structured charge vectors `v_r` such that

```text
g_0=||v_0||^2-2||v_1||^2+<v_0,v_2>.                  (P22.5zz14)
```

All charge indices in these displays are modulo `n`.

Remove one factor `B_a` of sign `epsilon` from an even-minus word and
call the remainder `R`.  For each feature polynomial of `R`, write

```text
B_gamma(z)=sum_r b_r(gamma)z^r,
b_(-r)(gamma)=epsilon b_r(gamma),
H_r(R)=sum_gamma w_gamma[
          b_r(gamma)^2-b_(r-1)(gamma)b_(r+1)(gamma)]. (P22.5zz15)
```

Then

```text
g_0(full word)=2g_a(R)=2[H_a(R)-H_(a+1)(R)].          (P22.5zz16)
```

Consequently odd-level orbit-ring `GKS2*` is equivalent to the structured
cyclic outer-Turan assertion

```text
H_a(R)>=H_(a+1)(R)                                   (COTM)
```

whenever every radius in `R` is at most `a<m` and the inversion parity
of `R` is `epsilon`.  It is enough to use `(COTM)` after removing a
factor of largest radius.

**Proof.**  Under the bijection `u=zw,v=z/w`, the first two identities in
`(P22.5zz11)` give the positive orthogonal feature expansion

```text
P(zw,z/w)=sum_gamma w_gamma A_gamma(z)A_gamma(w).
```

For even minus parity, `a_(-r)=a_r`.  Uniform cyclic averaging gives

```text
E A=a_0,
E C_1 A=2a_1,
E C_1^2 A=2a_0+2a_2.
```

The corner integral is one quarter of the average of

```text
[C_1(z)-C_1(w)]^2 A_gamma(z)A_gamma(w).
```

Expanding the square therefore gives the summand in `(P22.5zz13)`.
Taking the local feature vectors to be orthonormal and tensoring them gives

```text
<v_r,v_s>=sum_gamma w_gamma a_r(gamma)a_s(gamma),
```

which proves `(P22.5zz14)`.

It remains to telescope one complete local interval string.  Fix a feature
polynomial `B=sum_r b_rz^r` of the remainder.  When `epsilon=+1`, add

```text
2q(B)+sum_(t=1)^a q(C_tB);
```

when `epsilon=-1`, add

```text
sum_(t=1)^a q(S_tB),
```

where `q(A)=a_0^2-2a_1^2+a_0a_2`.  Direct expansion and cancellation of
the consecutive `t`-terms gives, in both cases,

```text
2[(b_a^2-b_(a-1)b_(a+1))
  -(b_(a+1)^2-b_a b_(a+2))].                        (P22.5zz17)
```

The endpoint convention is automatic modulo odd `n`; in particular, when
`a=m-1`, one has `b_(m+1)=b_(-m)=epsilon b_m`.  Summing
`(P22.5zz17)` over the remainder features proves the last equality in
`(P22.5zz16)`.  The first equality is the partial-character identity:
exchange parity of `R` is exactly the sign of the removed factor, so
adjoining it pairs the two equal contributions to the coefficient of
`B_a`.  Finally, removing a largest radius makes every remaining radius at
most `a`, proving the equivalence with `(COTM)`.  QED.

The form `(COTM)` is the odd finite, unit-spaced analogue of the ordinary
outer monotonicity target `(OTM)` in Proposition 15.  As there, the Gram
matrix alone is insufficient; the vectors come from tensor products of
complete nested interval strings.  The gain is that all affine wrapping is
now confined to the cyclic charge convention in `(P22.5zz15)`.  The
`turan-check` mode of `search_su2_odd_orbit_gks.cpp` independently expands
both sides of `(P22.5zz17)` as exact integer quadratic forms.  It passes all
9,900 radius/parity identities through rank one hundred, including every
cyclic endpoint case; see `certificates/su2_odd_orbit_gks.log`.

The first non-group orbit rank beyond Fibonacci can be settled completely.

**Proposition 22H3N (uniform `GKS2*` for the level-five orbit ring).**
The rank-three orbit ring `O_5` satisfies full `GKS2*` for signed words
of arbitrary length.  Consequently the corner `(SCW0)` holds at level
five for every plus word.  Together with the general partial-character
equivalence, this also proves the full boundary column `(SCB)` at level
five.

**Proof.**  Use the even-lift basis `1,B_1,B_2`.  Its fusion relations are

```text
B_1^2=1+B_1+B_2,       B_1B_2=B_1+B_2,
B_2^2=1+B_1.                                      (P22.5zz18)
```

Put `x=B_1` and `y=B_2=x^2-x-1`.  The three spectral values
`x_1>x_2>x_3` are the roots of

```text
f(X)=X^3-2X^2-X+1.
```

The cyclic model gives their masses

```text
w_i=(3-x_i)/7.                                    (P22.5zz19)
```

Direct rational sign checks give the isolating intervals

```text
2<x_1<9/4,       11/20<x_2<2/3,       -1<x_3<-4/5. (P22.5zz20)
```

Indeed `f(2)<0<f(9/4)`, `f(11/20)>0>f(2/3)`, and
`f(-1)<0<f(-4/5)`; the three disjoint intervals contain all three roots.

The map `T(X)=2+X-X^2` permutes these roots cyclically, since

```text
f(T(X))=f(X)(-X^3+X^2+2X-1),
```

and the intervals identify the images.
Since `y_i=1-T(x_i)`, after choosing the cyclic orientation one has

```text
y_1=1-x_3,       y_2=1-x_1,       y_3=1-x_2.       (P22.5zz21)
```

Set

```text
A=x_1-x_2,       B=x_2-x_3,       C=x_1-x_3=A+B,
alpha=A/B,       beta=C/A.
```

The root sum is two.  Since `x_2<2/3`, it follows that `A>B`, hence
`alpha>1` and `C/B>2`.  The sharper intervals in `(P22.5zz20)` give

```text
B/A>27/34>3/4,       beta>61/34>7/4,               (P22.5zz22)
w_1/w_2>1/4,         w_1/w_3>3/16,
w_2/w_3>7/12.                                       (P22.5zz23)
```

Apply the disjoint-support reduction.  There are only two nontrivial labels,
so an irreducible signed word has the form

```text
[B_1(x)+sigma B_1(y)]^p
[B_2(x)+tau B_2(y)]^q.                              (P22.5zz24)
```

All-plus words are positive by the fusion rule, and odd exchange parity
vanishes.  In every remaining even-parity case, diagonal spectral pairs
vanish and exchange of the two nodes gives

```text
J/2=sum_(i<j) w_iw_j
       (x_i+sigma x_j)^p(y_i+tau y_j)^q.             (P22.5zz24a)
```

If `sigma=-1,tau=+1`, take `p>=2` even.  Even `q` is
pointwise nonnegative on every unordered spectral pair.  For odd `q`, the
only negative pair is `{2,3}`; the positive pair `{1,3}` dominates it
because

```text
[w_1w_3 C^p x_1^q]/[w_2w_3 B^p |x_3|^q]
 =(w_1/w_2)(C/B)^p(x_1/|x_3|)^q>1.                  (P22.5zz25)
```

Here `w_1/w_2>1/4`, `C/B>2`, and `x_1/|x_3|>2`.

If `sigma=+1,tau=-1`, take `q>=2` even.  Odd `p` is the only
non-pointwise case.  The negative pair is again `{2,3}`, and now the
positive pair `{1,2}` dominates since

```text
[w_1w_2(2-x_3)^p C^q]/
 [w_2w_3(x_1-2)^p A^q]
 =(w_1/w_3)[(2-x_3)/(x_1-2)]^p(C/A)^q>1.             (P22.5zz26)
```

Indeed the first weight ratio is greater than `3/16`, while the middle
ratio is greater than eight.

It remains to put both labels in the minus support.  If `p,q` are even,
every pair is nonnegative.  If both are odd, write
`p=2r+1,q=2s+1` and put

```text
L_1=w_1w_2AC,       L_2=w_1w_3CB,       N=w_2w_3BA.
```

Orthogonality of the distinct characters `B_1,B_2` is exactly the
two-factor identity

```text
N=L_1+L_2.                                          (P22.5zz27)
```

After division by the negative pair's extra factor `B^(2r)A^(2s)`, the
two positive multipliers are

```text
M_1=alpha^(2r) beta^(2s),
M_2=alpha^(2(r-s)) beta^(2r).                        (P22.5zz28)
```

When `r>=s`, both are at least one and `(P22.5zz27)` finishes the
comparison.  If `s=r+d` with `d>=1`, discard the common factors at least
one.  It is enough to prove

```text
L_1 beta^(2d)+L_2 alpha^(-2d)>=L_1+L_2.
```

But

```text
(1-alpha^(-2d))/(beta^(2d)-1)
 <1/(beta^2-1)<16/33
 <7/12<(w_2/w_3)alpha=L_1/L_2,                       (P22.5zz29)
```

which is precisely that inequality.  Words containing only one minus label
to an even power are pointwise nonnegative.  These cases exhaust every
support-disjoint word.  Proposition 5 then restores words with overlapping
sign supports by induction on their length, proving full `GKS2*`.

Lemma 22H3J transfers the scalar statement to `(SCW0)`.  Proposition 21
applied to `O_5` gives all partial-character coefficients, and
`(P22.5zy)` transfers those to `(SCB)`.  QED.

Two independent exact C++ scans cross-check this all-length proof.  The
support-disjoint total-degree scan through degree five hundred tests 502,002
signed words and 1,505,005 partial-character coefficients.  The rectangular
rank-three scan tests every exponent pair `0<=p,q<=500`, all four sign
choices, and all three boundary targets: 1,004,004 words and 3,012,012
coefficients.  Both pass.  These finite scans validate the spectral formulas
and their implementation; Proposition 22H3N, not the bounds of the scans, is
the uniform argument.

The companion exact analyzer
`character_ring_iter/analyze_su2_odd_orbit_rank3_cone.cpp` also explains why
the preceding pair-dominance proof is useful.  In the standard six-dimensional
symmetric-square fusion basis the two plus operators have nonnegative
matrices, but each of the two squared-minus operators and the mixed-minus
operator already has negative entries.  Thus the ordinary symmetric-square
orthant is not a common invariant cone even at orbit rank three.  This is a
no-go statement for that naive cone, not an obstruction to another cone.
Full transcripts are in `certificates/su2_odd_orbit_gks.log`.

The cyclic Turan reduction has a closed one-step shell formula.  It does not
by itself prove `(COTM)`, but it identifies exactly which extra inequalities
an invariant cone must retain.

**Lemma 22H3O (cyclic shell update).**  In the notation of Lemma 22H3M put

```text
K_(r,s)(R)=sum_gamma w_gamma b_r(gamma)b_s(gamma),
W_(r,d)(R)=K_(r-d,r+d)-K_(r-d-1,r+d+1).              (P22.5zz30)
```

If `R'` is obtained from `R` by adjoining a factor of radius `d` and sign
`delta`, then

```text
H_r(R')=sum_(t=-d)^d H_(r+t)(R)+delta W_(r,d)(R),     (P22.5zz31)

H_r(R')-H_(r+1)(R')
 =H_(r-d)(R)-H_(r+d+1)(R)
   +delta[W_(r,d)(R)-W_(r+1,d)(R)].                   (P22.5zz32)
```

Consequently `(COTM)` for all words is equivalent to the two-sided shell
family

```text
H_(a-d)(R)-H_(a+d+1)(R)
 >=abs[W_(a,d)(R)-W_(a+1,d)(R)]                       (CSL)
```

for every `d<=a<m` and every remainder whose radii are at most `a`.

**Proof.**  For one local feature, multiplication by `C_t` replaces
`b_r` by `b_(r-t)+b_(r+t)`, while multiplication by `S_t` replaces it by
`b_(r-t)-b_(r+t)`.  In the update of `H_r`, the two diagonal terms give
`H_(r-t)+H_(r+t)`.  The two cross terms are

```text
W_(r,t)-W_(r,t-1).
```

They telescope from `t=1` to `d`.  The weight-two constant plus feature
supplies the remaining `H_r`, giving `(P22.5zz31)` for either sign.
Subtracting its translate by one telescopes the `H` sum and proves
`(P22.5zz32)`.  Adjoining the two possible signs gives the two inequalities
equivalent to `(CSL)`; conversely, remove any one factor from a nonempty
word and use `(P22.5zz32)`.  The empty word is immediate.  QED.

At orbit rank four, the spectral-pair problem is still finite enough for an
exact all-exponent transport certificate.  The elementary transport lemma
used by that certificate is recorded first.

**Lemma 22H3P (fractional exponential transport).**  Let

```text
F(r)=sum_(p in P)c_p product_a lambda_(p,a)^r_a
     -sum_(q in N)d_q product_a mu_(q,a)^r_a,          (P22.5zz33)
```

where all coefficients and bases are nonnegative and `r_a` range over the
nonnegative integers.  Join `q` to `p` when
`lambda_(p,a)>=mu_(q,a)` for every `a`.  If for every subset `X` of `N`,

```text
sum_(q in X)d_q<=sum_(p in Gamma(X))c_p,               (P22.5zz34)
```

then `F(r)>=0` for all `r`.

**Proof.**  The capacitated Hall theorem gives nonnegative transports
`theta_(q,p)`, supported on the displayed edges, whose row sums are `d_q`
and whose column sums are at most `c_p`.  On every edge the positive
monomial dominates the negative monomial coordinatewise for all `r`.
Multiply by `theta_(q,p)`, sum, and use the unused positive column capacity.
QED.

**Proposition 22H3Q (uniform `GKS2*` for the level-seven orbit ring).**
The rank-four orbit ring `O_7` satisfies full `GKS2*` for signed words of
arbitrary length.  Consequently `(SCW0)` and the full simple-current
boundary column `(SCB)` hold at level seven.

**Proof.**  Put `x_j=2cos(2 pi j/9)`, `1<=j<=4`.  The three values other
than `x_3=-1` are the roots, in decreasing order, of

```text
f(X)=X^3-3X+1.                                       (P22.5zz35)
```

In the even-lift basis `1,B_1,B_2,B_3`, set

```text
b_1(X)=1+X,
b_2(X)=X^2+X-1,
b_3(X)=X^3+X^2-2X-1,       w_j=(2-x_j)/9.             (P22.5zz36)
```

These are respectively the three orbit characters and their spectral
masses.  After the support-disjoint reduction, write a word as

```text
product_(a=1)^3[B_a(x)+s_a B_a(y)]^p_a,
s_a in {+1,-1},                                      (P22.5zz37)
```

omitting labels of exponent zero.  All-plus words are positive by fusion,
and odd exchange parity vanishes.  In every remaining even-parity chamber
the diagonal spectral terms vanish and

```text
J=2 sum_(i<j)w_iw_j
       product_(a=1)^3[b_a(x_i)+s_a b_a(x_j)]^p_a.    (P22.5zz38)
```

Fix the support, the signs, and the parity of each positive exponent.  Put
`p_a=p_a^0+2r_a`, where `p_a^0` is one or two and `r_a>=0`.  The signs of
the six summands in `(P22.5zz38)` are now fixed, while every increment of
`r_a` multiplies a summand by the nonnegative base

```text
lambda_((i,j),a)=[b_a(x_i)+s_a b_a(x_j)]^2.           (P22.5zz39)
```

Lemma 22H3P can therefore prove a whole exponent orthant from finitely many
base comparisons and Hall inequalities.

Here is the complete exact coverage.  There are nineteen even-parity
chambers containing a minus label in which every nonzero spectral pair is
already positive.  The remaining thirty chambers have negative pairs.  At
`r=(0,0,0)`, their thirty values are checked directly in the integral
fusion ring.  For every nonempty support
`E={a:r_a>0}`, first write `r_a=1+u_a` on `E`.  Of the resulting regimes,
161 satisfy Lemma 22H3P directly.

Four one-variable exceptions have signed minimum-exponent vectors

```text
(-1, 0,-1),       ( 0,+1,-2),
(-2,+1,-2),       (-1,+2,-1),                         (P22.5zz40)
```

and `E={3}`.  Peeling respectively two, one, one, and one values of `u_3`
leaves a tail satisfying Lemma 22H3P; the five peeled values are checked in
the exact fusion ring.

Six further regimes are covered by splitting the positive quadrant into
`u_1>=u_3` and `u_3>=u_1`.  In the first cone replace the two bases
`lambda_1,lambda_3` by `lambda_1 lambda_3,lambda_1`; in the second use
`lambda_1 lambda_3,lambda_3`.  Lemma 22H3P then applies.  Three last regimes
use the same two cones after common shifts four, three, and two.  The
complement consists of eighteen one-variable boundary strips.  Each strip
has a transported tail; its twenty-four peeled lattice points are checked
exactly.  These cases exhaust every support of `r`, so `(P22.5zz38)` is
nonnegative for arbitrary exponents.

For reproducibility, the strict verifier
`character_ring_iter/verify_su2_odd_orbit_rank4_transport.cpp` performs all
of the preceding tests with rational interval arithmetic.  It isolates the
three cubic nodes in

```text
1532088886237/10^12 < x_1 < 1532088886238/10^12,
 347296355333/10^12 < x_2 <  347296355334/10^12,
-1879385241572/10^12 < x_4 < -1879385241571/10^12,    (P22.5zz41)
```

checking the signs of `f` at both endpoints and the nonvanishing of `f'`.
Every dominance comparison and every one of the finitely many Hall
inequalities is then certified by ordered rational enclosures; no
floating-point decision is used.  The exact fusion calculation handles all
base and peeled points.  Its coverage line is

```text
chambers=30 pointwise_chambers=19 exact_base_points=30
direct_regimes=161 exact_peeled_points=5 tail_regimes=4
two_cone_regimes=6 shifted_two_cone_regimes=3
strip_regimes=18 strip_peeled_points=24
unresolved_regimes=0 result=PASS.
```

Thus full `GKS2*` holds in `O_7`; Proposition 5 restores overlapping sign
supports.  Lemma 22H3J transfers the scalar statement to `(SCW0)`.
Proposition 21 applied to `O_7`, followed by `(P22.5zy)`, gives the whole
column `(SCB)`.  QED.

This enlarges the all-length orbit theorem from ranks three to four, but it
is still a fixed-rank certificate.  The uniform odd-level problem is to
replace its finite spectral fan by a rank-independent transport, or
equivalently to prove the shell cone `(CSL)` directly.

The deterministic C++ diagnostic
`character_ring_iter/search_su2_odd_orbit_tropical.cpp` tests the most
immediate obstruction to such a transport.  For a sampled support, sign,
parity, and nonnegative exponent direction, it asks whether a negative
spectral pair is the unique maximizer of the logarithmic squared-base
functional.  Such a pair would dominate after scaling that direction and
would produce an eventual counterexample.  One million trials at each
orbit rank from five through fifteen, with direction coordinates at most
fifty, find no such pair.  This is bounded evidence only: it neither checks
every direction nor proves the required convex transport.  The transcript
is `certificates/su2_odd_orbit_tropical.log`.



The first rank-five transport exception has a uniform one-variable
description.  It converts that chamber into one sharp coefficient
log-concavity problem whose parameters no longer depend on a spectral fan.
The recursive transport analyzer still fails there after all primitive
two-dimensional rays with coordinates at most one hundred and every common
shift through one hundred.  Its transcript is
`certificates/su2_odd_orbit_transport_fan_rank5.log`.  This excludes that
pairwise Hall ansatz; it is not a negative orbit-ring value.

**Lemma 22H3R (top-orbit ray and separated-pole Turan reduction).**  At odd
level `k=2m-1`, let `Y=B_(m-1)=A_1` be the top even-lift orbit generator and
put `n=2m+1`.  For `M>=1` and `q>=0`, set

```text
J_(M,q)=< [Y(x)-Y(y)]^(2M)[Y(x)+Y(y)]^q >_(rho_m tensor rho_m).
```

Write `C(v)=v+v^(-1)` and

```text
A_(M,r)=(1/n)sum_(v^n=1)[4-C(v)^2]^M C(v)^r.        (P22.5zz42)
```

Then

```text
J_(M,q)=(1/2)[A_(M,q)A_(M,q+4)-A_(M,q+2)^2]  (q even),
J_(M,q)=(1/2)[A_(M,q+2)^2-A_(M,q)A_(M,q+4)]  (q odd). (P22.5zz43)
```

The first line is nonnegative pointwise.  For `q=2h+1`, define

```text
t_ell=4cos^2(ell pi/n),                    1<=ell<=m,
P_m(T)=U_(2m)(sqrt(T)/2)=product_ell(T-t_ell),
Q_m(z)=z^m P_m(1/z)=product_ell(1-t_ell z)
      =sum_(j=0)^m (-1)^j binom(2m-j,j)z^j,          (P22.5zz44)

F_(m,M)(z)=(4z-1)^(M-1)/Q_m(z)=sum_(N>=0)c_N z^N,
c_N=0                                                (N<0).
```

With `N=M+h-m`, one has the exact identity

```text
J_(M,2h+1)=2[c_(N+1)^2-c_N c_(N+2)].                (P22.5zz45)
```

Consequently full `GKS2*` on the one-generator top-orbit ray, uniformly in
the odd level, is equivalent to the separated-pole coefficient inequalities

```text
c_(N+1)^2>=c_N c_(N+2),
F_(m,M)(z)=(4z-1)^(M-1)/product_ell(1-t_ell z),
0<t_ell<4.                                          (SPT)
```

Only the displayed Chebyshev pole sets and indices `N=M+h-m` are required;
`(SPT)` is not being asserted here for arbitrary pole sets.  Its numerator
zero `1/4` lies strictly before every denominator pole `1/t_ell`.  Thus the
rank-five obstruction is a concrete sign-regularity problem for one rational
generating function, rather than a failure of the orbit inequality.

**Proof.**  Since the complete interval on `C_n` sums to zero away from the
identity,

```text
Y(u)=beta_(m-1)(u)=-(u^m+u^(-m)).
```

The automorphism `v=u^m`, whose inverse is `u=v^(-2)`, changes the orbit
measure into

```text
<f(Y)>_(rho_m)
 =(1/(2n))sum_(v^n=1)[4-C(v)^2]f(-C(v)).             (P22.5zz46)
```

For two roots put `v=ab,w=a/b`; this is a bijection because `n` is odd.
The elementary identities

```text
C(ab)+C(a/b)=C(a)C(b),
C(ab)-C(a/b)=(a-a^(-1))(b-b^(-1)),

[4-C(ab)^2][4-C(a/b)^2]=[C(a)^2-C(b)^2]^2
```

turn the double sum into a rank-two determinant.  Expanding its last square
gives `(P22.5zz43)`; the two factors
`(a-a^(-1))^(2M)` contribute the same sign, while the sign of the plus
factor is `(-1)^q`.

It remains to evaluate the odd cyclic moments.  Squaring is an automorphism
of `C_n`, its inverse exponent is `m+1`, and `C_(m+1)=C_m`.  Hence

```text
A_(M,2h+1)
 =(1/n)sum_(v^n=1) C_m(v)[2-C(v)]^M[2+C(v)]^h.
```

The identity root contributes zero.  On the other inverse pairs, ordered by
`ell=1,...,m`, one has

```text
C_m(v)=(-1)^ell sqrt(t_ell),
1/P_m'(t_ell)=(-1)^(ell+1)sqrt(t_ell)(4-t_ell)/n.
```

Therefore

```text
A_(M,2h+1)
 =-2 sum_(ell=1)^m
      t_ell^h(4-t_ell)^(M-1)/P_m'(t_ell)
 =-2 [t_1,...,t_m]{T^h(4-T)^(M-1)}.                 (P22.5zz47)
```

The last expression is an ordinary divided difference.  Expanding the
power of `4-T` and using

```text
[t_1,...,t_m]T^d=[z^(d-m+1)]1/Q_m(z)
```

with coefficients of negative degree set to zero gives

```text
[t_1,...,t_m]{T^h(4-T)^(M-1)}=c_(M+h-m).            (P22.5zz48)
```

Substitution in the odd line of `(P22.5zz43)` proves `(P22.5zz45)`.
The even line was already pointwise nonnegative, so the asserted equivalence
with `(SPT)` follows.  QED.

**Lemma 22H3S (separated-pole Toeplitz Turan theorem).**  Let `a>0`, let

```text
0<x_1<...<x_m<a,
R(z)=(az-1)^p/product_(i=1)^m(1-x_i z)=sum_(r>=0)c_r z^r,
```

where `p` is a nonnegative integer.  If

```text
r>=max(0,p+1-m),
```

then

```text
c_(r+1)^2-c_r c_(r+2)>=0.                           (P22.5zz49)
```

The inequality is strict for `m>=2`; for `m=1` it is an equality.  Thus
`(SPT)` holds for arbitrary distinct separated poles, not only for the
Chebyshev pole sets in Lemma 22H3R.

**Proof.**  The case `m=1` is the geometric tail following the polynomial
part, so assume `m>=2`.  Put

```text
q=r+m-1-p>=0,       f(X)=X^q(a-X)^p,
D_j=[x_1,...,x_m]{X^j f(X)}.
```

Partial fractions are valid beyond the polynomial part and give

```text
c_(r+j)=D_j,                                      j=0,1,2. (P22.5zz50)
```

Indeed, the coefficient of `(1-x_i z)^(-1)` is
`(a-x_i)^p x_i^(m-1-p)/P'(x_i)`, where
`P(X)=product_i(X-x_i)`.

Let `V=product_(i<j)(x_j-x_i)>0`, and let `E` be the determinant whose
rows, evaluated at the columns `x_1,...,x_m`, are

```text
f, Xf, 1, X, ..., X^(m-3).                          (P22.5zz51)
```

The last block is absent when `m=2`.  Expanding `(P22.5zz50)` in barycentric
form, or equivalently expanding `E` in its first two rows, gives

```text
D_1^2-D_0D_2=E/V.                                  (P22.5zz52)
```

For completeness, the contribution of a pair `i<j` to either side after
multiplication by `V` is obtained from

```text
-(x_i-x_j)^2/[P'(x_i)P'(x_j)]
 =1/product_(k notin {i,j})[(x_i-x_k)(x_j-x_k)].
```

It remains only to determine the sign of `E`.  Order the functions exactly
as in `(P22.5zz51)`.  Their first two leading Wronskians are `f` and `f^2`.
More generally, for `1<=s<=m-2`, expansion in the polynomial columns gives

```text
W(f,Xf,1,X,...,X^(s-1))
 =product_(j=0)^(s-1)j!
  {(s+1)[f^(s)]^2-s f^(s-1)f^(s+1)}.               (P22.5zz53)
```

Set `g=f^(s-1)`.  Since `f=X^q(a-X)^p` is real-rooted, so is every
derivative `g`, and Laguerre's elementary identity gives

```text
(g')^2-gg''>=0.
```

The braces in `(P22.5zz53)` are

```text
(g')^2+s[(g')^2-gg''].
```

They are strictly positive on `(0,a)`.  To see the only possible equality
case directly, all non-endpoint roots of every derivative of
`X^q(a-X)^p` are simple by strict Rolle interlacing.  Hence `g=g'=0`
cannot occur in the open interval, while at a point with `g'=0` and
`g!=0` the Laguerre expression is strict.  Also
`deg f=p+q=r+m-1>=m-1`, so none of the derivatives used above is constant.

Thus every leading Wronskian of `(P22.5zz51)` is positive on `(0,a)`.
The elementary Wronskian criterion (successive Rolle induction) says that
such an ordered family is an extended complete Chebyshev system.  Therefore
its evaluation determinant is positive whenever
`x_1<...<x_m`; that is, `E>0`.  Now `(P22.5zz52)` proves
`(P22.5zz49)`.  QED.

**Corollary 22H3T (full top-orbit ray).**  The one-generator top-orbit ray
of every odd fusion level satisfies `GKS2*` for arbitrary word length.

**Proof.**  In Lemma 22H3R take `a=4`, `p=M-1`, `r=N=M+h-m` and
`x_i=t_i`.  If `N>=0`, then

```text
N>=p+1-m,
```

so Lemma 22H3S proves `(SPT)`.  If `N<=-2`, both `c_N` and `c_(N+1)`
vanish; if `N=-1`, the Turan expression is `c_0^2`.  Hence the odd-plus
line of `(P22.5zz43)` is nonnegative for every `M,h`.  Its even-plus line
was already pointwise nonnegative.  QED.

The strict C++ analyzer `character_ring_iter/analyze_su2_tadpole_ray.cpp`
checks `(P22.5zz43)--(P22.5zz48)` by three independent integer calculations:
direct tadpole transfer, unfolding to the bipartite `A_(2m)` path, and
cyclic Laurent coefficients.  Its inequalities are still bounded checks;
the all-index conclusion is Corollary 22H3T, not an extrapolation from those
checks.  The reduction and Lemma 22H3S replace the first `O_9` spectral-fan
exception by a rank-independent Chebyshev-system argument.  The exact
transcript through rank twenty, twenty minus pairs, and plus exponent
thirty-nine is
`certificates/su2_tadpole_ray.log`.

The independent strict verifier
`character_ring_iter/verify_separated_pole_turan_ect.cpp` evaluates the
rational coefficients directly, computes the determinant `(P22.5zz51)` by
exact Bareiss elimination, and checks every leading-Wronskian identity
`(P22.5zz53)`.  Its bounded audit through separation nine, numerator power
ten, and coefficient degree sixteen comprises `77,518` coefficient checks,
`72,898` generalized-determinant identities, and `1,532,426` Wronskian
identities.  The transcript is
`certificates/separated_pole_turan_ect.log`.

Feeding Corollary 22H3T back into the rank-five chamber analyzer exposes
the next obstruction rather than merely renaming the old one.  At `O_9`
the former support-56 exception has odd exponents `p,r` on the two minus
labels `B_1,B_4`.  Since `Y=B_4` and

```text
B_1(x)-B_1(y)=[Y(x)-Y(y)][Y(x)+Y(y)],
```

that entire chamber is a top-orbit word with even difference exponent and
is proved by Corollary 22H3T.  The first still-untransported regime is

```text
[B_1(x)-B_1(y)]^(4+2u)[B_2(x)+B_2(y)]
[B_4(x)-B_4(y)]^(4+2v),                 u,v>=0.       (P22.5zz54)
```

Here `Y=B_4` and, on the top-orbit spectrum,

```text
B_1=Y^2-1,             B_2=Y^4-3Y^2+1=:Z.
```

Thus `(P22.5zz54)` is exactly the two-parameter partial-character target

```text
< [Y(x)-Y(y)]^(8+2u+2v)[Y(x)+Y(y)]^(4+2u)
   [Z(x)+Z(y)] > >=0.                                 (TP2)
```

It is twice the `B_2` partial coefficient of the even-plus top-orbit word.
The scalar theorem 22H3T does not by itself control that coefficient.  The
modified analyzer fails pairwise Hall transport for `(TP2)` after all
primitive rational two-dimensional rays with coordinates at most one
hundred and all common shifts through one hundred.  This is again a no-go
for that ansatz, not a negative value.  The transcript is
`certificates/su2_odd_orbit_transport_fan_rank5_after_spt.log`.

**Lemma 22H3U (three-feature partial top-orbit reduction).**  There is an
exact rank-independent reduction of this next target.  Put

```text
M=4+u+v,      q=4+2u,      a_j=A_(M,q+2j),   0<=j<=4,
b_j=a_(j+1)-2a_j,                              0<=j<=2,
e_j=a_(j+2)-4a_(j+1)+2a_j,                    0<=j<=2. (P22.5zz55)
```

Then the left side of `(TP2)` is

```text
(a_0a_2-a_1^2)
 +(1/2)(b_0b_2-b_1^2)
 +(1/2)(e_0e_2-e_1^2).                         (TPD) (P22.5zz56)
```

**Proof.**  After the root-of-unity change `v=ab,w=a/b` used in Lemma
22H3R,

```text
Z(-C(ab))+Z(-C(a/b))
 =2+C_2(a)C_2(b)+C_4(a)C_4(b).                 (P22.5zz57)
```

For any one of the three features `R=1,C_2,C_4`, multiplication by
`[C(a)^2-C(b)^2]^2` contributes one half of

```text
L_0(R)L_4(R)-L_2(R)^2,
L_j(R)=(1/n)sum_(z^n=1)(4-C(z)^2)^M C(z)^(q+j)R(C(z)).
```

The coefficient two of the constant feature and the identities
`C_2=C^2-2`, `C_4=C^4-4C^2+2` give `(P22.5zz56)`.

There is a shorter equivalent form.  Put

```text
R_1(T)=T-2,
R_2(T)=T^2-4T+2,
R_3(T)=T^3-6T^2+9T-2,
Lambda(G)=(1/n)sum_(z^n=1)
             (4-C(z)^2)^M C(z)^q G(C(z)^2).
```

The elementary Christoffel--Darboux identity

```text
2+R_1(T)R_1(U)+R_2(T)R_2(U)
 =[R_3(T)R_2(U)-R_2(T)R_3(U)]/(T-U)                (P22.5zz58)
```

turns `(TPD)` into the single cross determinant

```text
2(TPD)=Lambda(TR_3)Lambda(R_2)
       -Lambda(R_3)Lambda(TR_2).                    (P22.5zz59)
```

To obtain it, write every ordinary `2 by 2` moment determinant as one half
of its double-sum form, insert `(P22.5zz58)`, and exchange `T,U` in the
second term.  This also proves `(P22.5zz59)`.
QED.

Unlike Lemma 22H3S, the last two determinants are moments of signed
features, so their separate positivity is false in general; `(TPD)` is the
new combined target.  The tadpole analyzer independently compares `(TPD)`
with direct orbit transfer in `7,056` exact cases through rank twenty and
`0<=u,v<=20`, all nonnegative.  This is bounded evidence, not yet an
all-index proof.

**Lemma 22H3V (stable and no-alias `TP2`).**  The target `(TP2)` is
strictly positive in the continuous cyclic limit.  At finite orbit rank
`m`, it is also positive whenever

```text
M+r+4<2m+1,
M=4+u+v,       r=2+u.                               (P22.5zz60)
```

In particular, every fixed pair `u,v` is proved for all sufficiently large
odd levels; the base pair `u=v=0` is already in the no-alias range at
`O_9`.

**Proof.**  Let `Lambda_infinity` be the uniform circle limit of `Lambda`.
Under `T=2+2cos(theta)`, its tilted measure is, up to a positive constant,

```text
T^(r-1/2)(4-T)^(M-1/2)dT,             0<T<4.
```

Put

```text
d=M-r,       s=M+r+1,
ell_j=Lambda_infinity(R_j)/Lambda_infinity(R_0).
```

The beta Pearson identity and the Chebyshev derivative identity are

```text
Lambda_infinity(T(4-T)F')
 =Lambda_infinity([sT-(4r+2)]F),

T(4-T)R_j'=j(R_(j-1)-R_(j+1)).                      (P22.5zz61)
```

Together with `TR_j=R_(j+1)+2R_j+R_(j-1)`, they give

```text
(s+j)ell_(j+1)+2d ell_j+(s-j)ell_(j-1)=0,
ell_0=1,       ell_1=-d/s.                           (P22.5zz62)
```

Set

```text
A=2d^2-s^2+s,       B=3s^2-3s-2-4d^2.
```

The first two needed coefficients are

```text
ell_2=A/[s(s+1)],
ell_3=dB/[s(s+1)(s+2)].                              (P22.5zz63)
```

If `J` denotes the right side of `(P22.5zz59)` divided by
`Lambda_infinity(R_0)^2`, two further uses of `(P22.5zz62)` give

```text
J=N/[s^2(s+1)^2(s+2)^2(s+3)(s-2)],                  (P22.5zz64)

N=6(s-2)(s+2)^2 A^2+10d^2(s+2)AB+4(s+3)d^2B^2.
```

The middle summand in this last presentation is not manifestly positive.
The required parameter substitution is decisive: `r=2+u`, `d=2+v`, hence
`s=7+2u+v`.  Direct expansion gives

```text
N=sum_(i=0)^7 u^i P_i(v),                            (P22.5zz65)

P_0=3353400+2915280v+985070v^2+170330v^3
    +18870v^4+1750v^5+100v^6,
P_1=7477920+5687192v+1623228v^2+220152v^3
    +16628v^4+1000v^5+40v^6,
P_2=7027704+4555912v+1061000v^2+107048v^3
    +4792v^4+120v^5,
P_3=3605232+1915712v+343408v^2+23296v^3+464v^4,
P_4=1089600+445248v+54912v^2+1920v^3,
P_5=193920+54144v+3456v^2,
P_6=18816+2688v,
P_7=768.
```

Every coefficient is strictly positive, while the denominator in
`(P22.5zz64)` is positive because `s>=7`.  Equations `(P22.5zz59)` and
`(P22.5zz64)` therefore prove the continuous assertion.

For the finite statement, write the circle variable after the squaring
automorphism.  The Laurent degree of
`T^r(4-T)^M R_j(T)` is at most `M+r+j`.  Every functional in
`(P22.5zz59)` has `j<=4`.  If `(P22.5zz60)` holds, averaging over the
`(2m+1)`st roots of unity extracts only the ordinary constant term, exactly
as continuous circle averaging does.  Thus the same positive determinant
applies.  For fixed `u,v`, condition `(P22.5zz60)` eventually holds as
`m` tends to infinity.  QED.

The strict C++ program `character_ring_iter/derive_continuous_tp2.cpp`
performs the expansion `(P22.5zz65)` independently in exact bivariate
`cpp_int` arithmetic.  Its transcript, including all thirty-five positive
coefficients, is `certificates/continuous_tp2.log`.  Lemma 22H3V leaves
only the finite cyclic-alias range

```text
M+r+4>=2m+1                                           (TP2-alias)
```

for this chamber.

The exact diagnostic `character_ring_iter/probe_cyclic_christoffel_minors.cpp`
tests the broader cross-minor family

```text
Lambda(TR_(j+1))Lambda(R_j)-Lambda(R_(j+1))Lambda(TR_j).
```

Through rank twenty and endpoint powers thirty, all `7,182` cases with
`j=2`, `r>=2`, and `M-r>=2` pass, including the alias range.  The stronger
`j=2` scan through rank twenty-five and both powers forty finds no negative
minor in any of its `80,688` cases, even outside the target parameter cone.
The still broader unrestricted pattern is false: the `j<=6` run finds
`7,664` negative minors
at higher indices or outside the target parameter cone.  Thus the evidence
supports exactly the needed theorem and rules out a blanket
connection-minor positivity claim.  The transcript is
`certificates/cyclic_christoffel_minors.log`.

Two larger standard packet cones still do not isolate that differential.
For the level-three word `Q=[2,2,2]`, the full affine packet is

```text
5 Psi_(0,0)-6 Psi_(1,1)+2 Psi_(2,0)+5 Psi_(2,2).      (P22.5zq)
```

Its constant-`r+s=2` inward prefix has value `2-6=-4`, although its row
coefficient is `2`.  Thus diagonal-prefix positivity is false.  Adding the
lower-degree reservoir gives the canonical two-layer packet
`Psi_(2,0)+Psi_(1,1)+Psi_(0,0)`, which has value one here, but that larger
cone also eventually fails.  The first bounded obstruction is

```text
k=5, Q=[1,1,3,3,4,4], (a,b)=(3,3):
two-layer packet=-2.                                  (P22.5zr)
```

The strict mode `finite-c2-wall-packets` finds this after checking every
smaller level/word in its bounds; all its `b=1` packets, which are exactly
the row coefficients `(P22.5zp)`, remain nonnegative before the obstruction.
Hence the required invariant cone is genuinely row-adapted: it is weaker
than character positivity, diagonal-prefix positivity, and the ordinary
two-layer packet inequalities.

The nearest-neighbor hypothesis is substantive.  A general label `q` gives
a fusion jump which may cross `x+y=k` without landing on it; reflecting only
the suffix would then split one irreducible `q`-block between the two
coordinates.  Thus Lemma 22H3F does not by itself extend to arbitrary plus
labels, but it supplies a uniform finite-level base cone on which a
block-level recoupling proof of `(SCBI)` can build.

The two reservoirs on the right of `(SCBI)` must be coupled, even in the
fully support-disjoint chamber and when the next label is strictly larger
than every prefix label.  The first exact example is

```text
k=5, Q=[2,2,2,2,2,3,3], q=4, a=1:
boundary interval sum=101,       lower term=249,
interior sink=327,               combined slack=23.    (P22.5zi)
```

Here the final minus terminals are `1,5`, the sink terminals are `4,4`, and
none occurs in `Q`; the appended `q=4` is also new.  Thus both `101<327` and
`249<327`, whereas their sum pays the sink.  The strict C++ mode
`search_su2_level_partial_character_cone simple-current-boundary` verifies
the exact recurrence and `(SCB)` in 296,517 ordered prefix/next/target
entries through level ten and seven prefix factors.  It finds 7,483
fully-disjoint strict-next instances in which both individual reservoirs
fail, while the combined recurrence always passes in the bounded domain.
The output is in `certificates/su2_simple_current_boundary.log`.  Hence a
proof of `(SCBI)` must transfer capacity between its two displayed positive
terms; separate domination is false.

Nor can that transfer be confined to one allocation of the plus factors and
its coordinate-reversed complement.  Decompose `C_Q` linearly into matrices

```text
x tensor y+y tensor x,
```

where `x` and `y` are the two subset fusion-multiplicity vectors.  At level
three, with

```text
Q=[1,1,1,1],       allocation mask=1,       q=2,       a=3,
```

the corresponding local `(SCBI)` terms are

```text
boundary=-1,       lower=-1,       sink=0,       slack=-2.  (P22.5zn)
```

Even summing all complementary allocation pairs of a fixed cardinality does
not repair the inequality: cardinality shell one in the same example has
slack `-8`.  A separate boundary failure occurs for
`Q=[1,1,2,2]`, shell two, with value `-2`.  The global sum is nevertheless
nonnegative.  Thus the required block switch must move between distinct
complement pairs and change the number of whole factors assigned to a
coordinate, just as the minuscule suffix reflection does.  The exact mode
`simple-current-boundary-local` records these obstructions in
`certificates/su2_finite_c2_boundary.log`.

For `Q=[2,2,2,2,3]` and `p=2`, exact fusion-path enumeration gives

```text
A_(1,2)=44,              A_(0,3)=43,              A_(0,1)=28. (P22.5y)
```

Thus the upper-branch-only strengthening `A_(1,p)<=A_(0,p+1)` already fails
in the first hard Plucker case, by one dimension.  The lower `p-1` branch
has ample supply and the combined coefficient has value `27`, but it must
be coupled to the merged row-syzygy space rather than omitted.  The strict
C++ enumerator `character_ring_iter/search_su2_fundamental_transfer.cpp`
finds 383 such upper-only failures in 1,261 words through nine plus factors,
label seven, and total degree twenty-eight, while `(P22.5x)` passes
throughout.
The exact output is in
`certificates/su2_fundamental_transfer.log`.  This is bounded evidence for
`(P22.5x)` and an exact obstruction to a one-branch proof.

A simultaneous sign gauge of the exterior fusion matrices does not prove
this target, even after the support-disjoint reduction removes plus label
one.  Write `s_(u,v)` for a prospective sign of the wedge state
`e_u wedge e_v`, with `u>v`.  Nonnegativity of `dGamma(N_2)` would force

```text
s_(3,0)=s_(1,0),              s_(3,2)=s_(3,0),
```

from its positive transitions `(1,0)->(3,0)` and `(3,0)->(3,2)`.
Nonnegativity of `dGamma(N_3)` would force

```text
s_(2,0)=s_(1,0),              s_(3,2)=-s_(2,0),
```

because the transitions `(2,0)->(1,0)` and `(2,0)->(3,2)` have signs
respectively positive and negative.  These four requirements contradict one
another.  The strict C++ signed-union-find diagnostic
`character_ring_iter/analyze_su2_additive_compound_sign.cpp` reproduces the
conflict both on the ordinary half-line and already on the finite level-three
fusion graph.  Hence Target 22H3 must use the ordered product or the complete
channel grouping; entrywise positivity of all gauged one-step compounds is
impossible.  The exact outputs and four-edge sign certificate are recorded
in `certificates/su2_additive_compound_sign.log`.

The same construction makes sense in every semisimple `SU(2)_k` fusion
category, with fusion-channel Hom spaces and the associator replacing
ordinary tensor reordering.  Thus a uniform complex whose degree-zero map is
`mu_-` and whose next chain group is assembled from the nonempty positive
channels would prove `(PCP_k)` directly.  Proposition `(7.2)` is exactly the
seven-factor, target-zero instance of this kernel bound.  The remaining
problem is no longer to guess the Euler characteristic: it is to construct
the Plucker/Temperley--Lieb differential and prove that it covers
`ker(mu_-)` without using more positive-channel dimension than `P_p`.

There is an exact tableau refinement of this target which does not require a
choice of Plucker straightening coefficients.  Adjoin a distinguished target
vertex `star` of degree `p` and choose any total order on the remainder
vertices together with `star`.  For a degree vector `d=(d_v)` of even total
degree, put `2D=sum_v d_v` and let `B(d)` be the set of integer vectors
`alpha=(alpha_v)` satisfying

```text
0<=alpha_v<=d_v,              sum_v alpha_v=D,
sum_(v<=w)(d_v-alpha_v)<=sum_(v<w)alpha_v   for every vertex w.   (P22.5)
```

The vector `alpha` records the multiplicity of each vertex in the top row of
a semistandard tableau of rectangular shape `(D,D)` and content `d`; the
bottom-row multiplicity is `d_v-alpha_v`.  The last inequalities say exactly
that, after both rows are sorted, every column is strictly increasing.
For odd total degree set `B(d)=emptyset`.

For a cut `T` not containing `star`, write `d_T` for the degree vector
obtained by zeroing all vertices outside `T`, and write
`d_(T^c star)` for its complementary degree vector, including `star`.  If
`alpha in B(d)`, call `T` admissible for `alpha` when

```text
alpha|_T in B(d_T),
alpha|_(T^c union {star}) in B(d_(T^c star)),             (P22.6)
```

and denote the family of such cuts by `F_alpha`.

**Proposition 22A (diagonal-tableau fibre identity).**  For every signed
remainder, target `p`, and total order of the vertices,

```text
dim U_(p,T)=number of alpha in B(d) for which T in F_alpha,       (P22.7)

g_p=sum_(alpha in B(d))
       sum_(T in F_alpha)(-1)^(|T intersect M|).                  (P22.8)
```

**Proof.**  The standard monomial basis of the multigraded bracket ring is
the set of semistandard two-row tableaux.  Once the content and the top-row
multiplicities are fixed, both sorted rows, and hence the tableau, are
unique.  This identifies the basis with `(P22.5)`.  It applies to
`I_T=Inv(tensor_(i in T)V_i)` and, after adjoining the target vertex, to
`C_(p,T)=Inv(V_p tensor tensor_(i notin T)V_i)`.

Take one tableau from each factor.  Their vertex supports are disjoint, so
adding their top-row multiplicity vectors gives an element `alpha` of
`B(d)`: every inequality in `(P22.5)` is the sum of the corresponding two
inequalities.  Conversely, for fixed `alpha` and `T`, the two factor vectors
can only be the restrictions in `(P22.6)`.  Thus a pair of factor tableaux
exists if and only if `T in F_alpha`, and it is then unique.  Summing over
`alpha` proves `(P22.7)`.  Substitution in `(P22.4)` proves `(P22.8)`.  QED.

The fibre identity suggests the following sufficient strengthening.  Suppose
that equal nontrivial remainder labels all carry the same sign, as allowed by
Proposition 5, and seek a total order for which every `alpha in B(d)` obeys

```text
number of T in F_alpha with |T intersect M| even
 >= number of T in F_alpha with |T intersect M| odd.             (P22.9)
```

This strengthening is false, even after optimizing the vertex order.

**Proposition 22B (order-independent fibrewise obstruction).**  Take

```text
remainder labels=[1,2,3,3,3,3],       target=3,
minus labels=[1,2],                    plus labels=[3,3,3,3].     (P22.10)
```

The remainder is support-disjoint.  For every total order of its six
vertices and `star`, some diagonal-tableau fibre violates `(P22.9)`.

**Proof.**  There are `7!=5,040` orders before quotienting by the four
identical plus vertices.  The strict integer enumerator constructs `B(d)`
from `(P22.5)`, tests both restrictions in `(P22.6)`, and exhausts every
order.  No order passes all fibres.  In an order maximizing the number of
odd fibres hit, the odd source dimension is `20`, the number of hit fibres
is `15`, and the nonempty even supply is `80`; one remaining bad fibre has
four odd and three even cuts.  The replay is
`certificates/su2_plucker_fibre_obstruction.log`.  QED.

Thus compensation must move between different tableau fibres.  For a fixed
order put

```text
o_alpha=number of odd cuts in F_alpha,
e_alpha=number of nonempty even cuts in F_alpha,
r_-=number of alpha with o_alpha>0.                              (P22.11)
```

In the diagonal Gelfand--Tsetlin degeneration, all odd basis products in one
fibre have the same monomial image and products in distinct fibres have
distinct images.  Its negative multiplication kernel therefore has
dimension

```text
kappa_GT=sum_alpha(o_alpha-indicator_(o_alpha>0)),
P_p=sum_alpha e_alpha.                                           (P22.12)
```

There is a concrete capacitated version of the required cross-fibre
compensation.  For a fixed order, let `O_alpha` be the odd cuts in
`F_alpha` and let `E_alpha` be its nonempty even cuts.  Whenever
`O_alpha` is nonempty, choose an anchor `A_alpha in O_alpha`.  Pair any

```text
l_alpha=min(|O_alpha|-1,|E_alpha|)
```

of the remaining odd cuts with distinct members of `E_alpha`; these are
the local payments.  Denote the unpaired odd cuts by `R_alpha`.  For every
nonempty even cut `S`, put

```text
u(S)=number of locally used pairs (alpha,S),
d(S)=number of (alpha,T) with T in R_alpha and A_alpha triangle T=S,
c(S)=dim U_(p,S).                                           (P22.13)
```

Here `triangle` denotes symmetric difference.  Since both `A_alpha` and
`T` are odd, `S` is even; since they are distinct, it is nonempty.

**Lemma 22C1 (local-plus-XOR capacity implication).**  If the order,
anchors, and local pairings can be chosen so that

```text
d(S)+u(S)<=c(S) for every nonempty even S,                   (HXC)
```

then `(AGT)` holds.

**Proof.**  The locally paired objects are distinct positive source-basis
objects `(alpha,S)`.  For fixed `S`, the remaining positive source contains
exactly `c(S)-u(S)` unused objects with cut `S`, irrespective of their
tableau fibres.  Condition `(HXC)` injects the `d(S)` residual relations
into those objects.  Therefore

```text
kappa_GT
 =sum_alpha(|O_alpha|-1)
 =sum_alpha l_alpha+sum_S d(S)
 <=sum_S u(S)+sum_S(c(S)-u(S))=P_p.
```

This is `(AGT)`.  QED.

**Target 22C2 (adaptive hybrid-capacity lemma).**  Prove that every
support-disjoint signed remainder and target admit a vertex order, anchors,
and local pairings satisfying `(HXC)`.  This target is strictly explicit:
it is a capacitated Hall matching between the special-fibre kernel
generators and the positive source-basis objects, with edges either inside
one fibre or in the symmetric-difference channel.

**Target 22C (cross-fibre Hall matching).**  Prove, for every
support-disjoint remainder and target, that some vertex order satisfies

```text
There is an order with kappa_GT<=P_p.                            (AGT)
```

Unlike `(P22.9)`, `(AGT)` holds in `(P22.10)`: the displayed best order has
`kappa_GT=20-15=5<=80`.

**Proposition 22D (aggregate degeneration transfer).**  Target 22C implies
`(UKS)`, hence `(PCP)` and the full ordinary `SU(2)` central-character `Q3`
theorem.

**Proof.**  The diagonal toric algebra is a flat degeneration of the bracket
ring, and the direct sum of odd multiplication maps specializes with it.
Matrix rank can only drop on the special fibre.  Hence the original negative
kernel dimension is at most `kappa_GT`, while every multigraded channel
dimension in `P_p` is constant in the flat family.  Target 22C gives
`kappa_p<=kappa_GT<=P_p`, which is `(UKS)`.  Proposition 22 and Proposition
19 finish the implication.  QED.

Target 22C is the current cumulative Plucker target; a fibre-preserving
involution is no longer viable.

Several natural strengthenings of Target 22C2 are false.

**Proposition 22E (exchange and fixed-fibre obstructions).**

1.  The sets `B(d)` are not, in general, bases of discrete polymatroids.
    For

```text
d=[1,1,2,1,1]
```

    they contain

```text
u=[1,1,0,1,0],                 v=[1,0,2,0,0].
```

    The only coordinate toward which the excess of `u` in coordinate two
    can move is coordinate three, but
    `u-epsilon_2+epsilon_3=[1,0,1,1,0]` is not in `B(d)`.  Thus neither the
    lattice-path-polymatroid exchange theorem nor its sorting Gröbner basis
    applies directly to `(P22.5)`.

2.  Pairwise XOR capacity can vanish even for two support-disjoint odd cuts
    in one fibre.  In the order with degree vector

```text
[1,1,1,3,2,2_star]
```

    take only the label-two remainder vertex negative.  The fibre

```text
alpha=[1,1,1,2,0,0]
```

    admits the two odd cuts with masks `19` and `25`, but their symmetric
    difference has mask `10`, carrying labels `[1,3]`; hence its invariant
    multiplicity and `c(10)` are zero.

3.  Even after restricting to the genuine source-deficit sector, no order
    need balance every fibre.  Take

```text
remainder=[1^-,2^+,2^+,2^+,2^+,3^+,4^+],       target=2.   (P22.14)
```

    Here `dim O_p=180` and `P_p=156`, so a rank contribution of at least
    `24` is genuinely needed.  Exhaustion of all `8!=40,320` vertex orders
    finds no order satisfying `(P22.9)`.  Nevertheless the order

```text
[1^-,2^+,3^+,4^+,2^+,2^+,2^+,2_star]
```

    has `r_-=74`, hence `kappa_GT=106<=156`, and its deterministic
    local-plus-XOR payment satisfies `(HXC)`.

**Proof.**  The first claim is the displayed failed base exchange.  For the
second, direct substitution in `(P22.5)--(P22.6)` gives the listed fibre and
cuts; `Inv(V_1 tensor V_3)=0`.  For the third, the strict integer enumerator
constructs every fibre for every order.  Its `gt-case` mode exhausts all
`40,320` permutations, while its hybrid mode counts every local use and
residual demand by cut mask and checks `(HXC)`.  The exact outputs are in
`certificates/su2_plucker_adaptive_capacity.log`.  QED.

The distinguished target can be removed from the combinatorics.  Give every
remainder vertex its prescribed sign `s_i` and give the target vertex the
exchange sign

```text
s_star=product_i s_i.                                  (P22.15)
```

The product of all vertex signs is then `+1`.

**Lemma 22F0 (root-free even-sign form).**  Regard a channel as an unordered
bipartition of the enlarged vertex set into two invariant-admissible sides.
Its sign is the product of the vertex signs on either side.  This sign is
well-defined, and the coefficient `g_p`, the dimensions `dim O_p,P_p`, and
the GT rank `r(prec)` are unchanged if any vertex is chosen as the target.

**Proof.**  Since the product of all signs is one, the two side products are
equal.  For a chosen target, exactly one side of every unordered bipartition
omits it; this is the cut `T` in `(P22.1)`.  The two invariant multiplicities
in `U_(p,T)` are the multiplicities of the two sides and their product is
symmetric under interchanging them.  Hence changing the target merely changes
which side represents the same unordered bipartition.  The empty channel is
the unique trivial bipartition and remains empty for every target.  Thus the
signed coefficient, the odd and nontrivial even source dimensions, and the
sets of tableau products in every ordered fibre are all unchanged.  QED.

**Corollary 22F0B (signed-class quotient).**  Vertices with the same pair
`(degree,sign)` are interchangeable, including the target after `(P22.15)`.
Consequently `r(prec)` is a function of the word of signed degree classes.
If class `c` has multiplicity `m_c`, every class word represents exactly
`product_c m_c!` physical orders.  Targets 22F1 and 22F2 may therefore be
checked on the multiset-permutation graph without changing either statement.

**Proof.**  Relabelling equal signed degrees preserves every box, queue,
cut parity, and tableau product.  If the relabelling moves the chosen target,
first apply Lemma 22F0.  The fibre size of the map from physical orders to
class words is the displayed constant, and adjacent unequal class swaps are
exactly the nontrivial edges of the quotient order graph.  QED.

There is an exact queue description of the order dependence.  Fix an order
and write

```text
beta_i=d_i-alpha_i.
```

For a cut `T`, colour its vertices one and the complementary vertices zero;
the target vertex has colour zero.  For `c=0,1` put

```text
H_c(j)=sum_(i<j, colour(i)=c) alpha_i
       -sum_(i<=j, colour(i)=c) beta_i.               (P22.16)
```

Thus a vertex first consumes its `beta_i` units from its colour queue and
then produces its `alpha_i` units.

**Lemma 22F (two-queue and adjacent-swap normal form).**

1.  A cut `T` belongs to `F_alpha` if and only if both colour queues in
    `(P22.16)` are nonnegative at every vertex and each returns to zero.

2.  Suppose positions `i,i+1` are interchanged.  Set

```text
h=sum_(j<i)(alpha_j-beta_j),
z=max(0,beta_i+beta_(i+1)-h).                          (P22.17)
```

    The Bender--Knuth bijection from the old `B(d)` to the `B(d')` for the
    interchanged degree sequence is

```text
alpha'_i     =alpha_(i+1)+z,
alpha'_(i+1) =alpha_i-z,                               (P22.18)
```

    with every other coordinate fixed.  It is an involution after the two
    degrees are interchanged.

**Proof.**  The first assertion is `(P22.5)` applied separately to the two
restrictions in `(P22.6)`: its prefix inequality says exactly that the
current queue can pay the current `beta_i`, and its total-degree condition
says that the queue returns to zero.

For the second, global validity at positions `i` and `i+1` gives

```text
h>=beta_i,
h-beta_i+alpha_i>=beta_(i+1).
```

Consequently `0<=z<=beta_(i+1)` and `0<=z<=alpha_i`.  After `(P22.18)` the
two new bottom multiplicities are

```text
beta'_i=beta_(i+1)-z,          beta'_(i+1)=beta_i+z.   (P22.19)
```

If `z=0`, then `h>=beta_i+beta_(i+1)`; if `z>0`, then
`h=beta_i+beta_(i+1)-z`.  In either case the two new prefix inequalities
follow immediately from `(P22.18)--(P22.19)`.  The combined top and bottom
totals are unchanged.  Applying the same formula to the interchanged pair
again gives the same `z` and recovers `alpha`.  This proves both validity
and bijectivity.  QED.

**Corollary 22F0A (fixed-universe rank).**  Let

```text
A(d)={alpha: 0<=alpha_i<=d_i and
             sum_i alpha_i=(sum_i d_i)/2}.             (P22.20)
```

For an order `prec`, `r(prec)` is the number of `alpha in A(d)` for which
there is an odd two-colouring whose two induced job sequences both satisfy
the queue condition `(P22.16)`.  In particular, adjacent orders can be
compared inside the same box slice `A(d)`; no transport of `alpha` is needed.

**Proof.**  If both colour queues are feasible and return to zero, their sum
is the feasible global queue, so `alpha in B_prec(d)`.  Lemma 22F then says
that the colouring is precisely an odd admissible cut.  The converse is the
same lemma.  QED.

For an order `prec`, let

```text
r(prec)=number of alpha in B_prec(d) admitting at least one odd cut,
R=max(0,dim O_p-P_p).                                  (P22.21)
```

Both dimensions in `R` are independent of the order, and `(AGT)` is exactly
`r(prec)>=R`.

**Target 22F1 (rank-threshold plateau lemma).**  In the adjacent-transposition
graph of vertex orders, every connected component of a level set
`r(prec)=r<R` has an adjacent edge to an order of strictly larger rank.

Target 22F1 implies Target 22C: move inside the current rank plateau to such
an edge and cross it.  Rank increases at most `|B(d)|` times, so the process
terminates at an order with `r>=R`.  This target is weaker than monotone
reachability of an `(HXC)` certificate.  The latter is already false for the
deterministic anchor/payment rule: for

```text
remainder=[1^+,1^+,2^+,2^+,3^+,4^-],       target=1,  (P22.22)
```

there are 32 trapped orders in the `7!`-order landscape.  This does not
obstruct `(AGT)`: here `dim O_p=25`, `P_p=18`, so `R=7`, while the best order
has rank `15` and kernel `10<=18`.  Every subthreshold order in this
landscape does satisfy Target 22F1.

The formulas `(P22.17)--(P22.18)` make Target 22F1 a local statement about two coloured
queues.  A proof must show that a subthreshold plateau cannot be closed under
all adjacent queue transfers; no choice of anchors or local payments occurs
in this reduced statement.

There is an even weaker averaged target.  If `o_alpha(prec)` is the number of
odd admissible cuts in one ordered fibre, then for every order

```text
dim O_p-r(prec)
 =sum_alpha[o_alpha(prec)-indicator_(o_alpha(prec)>0)]. (P22.23)
```

**Target 22F2 (average-order relation bound).**  With `n` remainder vertices,

```text
sum_prec sum_alpha
 [o_alpha(prec)-indicator_(o_alpha(prec)>0)]
 <=(n+1)! P_p.                                         (P22.24)
```

Equivalently,

```text
(1/(n+1)!) sum_prec r(prec)>=dim O_p-P_p.              (P22.25)
```

Target 22F2 implies Target 22C by averaging.  Unlike Target 22F1, it does not
require a monotone path or even a preferred order.  The strict `gt-average`
mode passed all `2,422` support-disjoint cases with at most seven factors,
labels and targets at most four, and total degree at most sixteen.  After the
signed-class quotient of Corollary 22F0B, both this mode and the independent
rank-threshold plateau mode passed all `3,201` cases at the same label and
factor bounds with total degree at most eighteen.  The exact replays and the
three obstruction outputs below are recorded in
`certificates/su2_gt_rank_average.log`.

Two tempting simplifications of `(P22.23)` are false.  First, the transports
`(P22.18)` do not define a symmetric-group action: for

```text
d=[1,1,1,1],                 alpha=[1,0,1,0],          (P22.26)
```

the two braid words at the first three positions give respectively
`[1,1,0,0]` and `[1,0,1,0]`.  Thus an averaged proof must work on the order
graph itself, not on path-independent Bender--Knuth orbits.  Second, replacing
the minimal relation count `o_alpha-1` by all odd pairs
`binomial(o_alpha,2)` is too expensive.  For

```text
remainder=[1^+,1^+,1^+,1^+,2^-],          target=0,   (P22.27)
```

one has `dim O_0=9`, `P_0=6`, and the sum of the ranks over all `6!=720`
orders is `2,160`, exactly the required average rank three.  But the sum of
the odd-pair collision counts is `7,488`, exceeding the available averaged
capacity `720*6=4,320`.  A proof of `(P22.23)` must therefore select a spanning
set of fibre relations rather than charge every colliding pair.

Nor can `(P22.24)` be proved separately for each raw `alpha` in the fixed
box slice of Corollary 22F0A.  The first exact obstruction is

```text
remainder=[1^-,1^-,2^+,3^+],             target=1,
alpha=[0,0,0,3,1].                                  (P22.28)
```

Across all `5!=120` physical orders this `alpha` carries fourteen excess odd
relations and no positive object with the same raw top vector.  The complete
case still satisfies `(AGT)`: `dim O_1=2`, `P_1=1`, and rank one suffices.
Thus an averaged proof must transfer capacity between different `alpha` as
well as between different orders.  Formula `(P22.18)` supplies exactly such
an adjacent transfer, even though `(P22.26)` shows it is path-dependent.

The path dependence can be retained rather than quotiented out.  Define the
Bender--Knuth state graph `K(d)` to have vertices

```text
(prec,alpha),                  alpha in B_prec(d),     (P22.29)
```

and join two vertices when their orders differ by one adjacent transposition
and their top vectors are related by `(P22.18)`.  Lemma 22F makes every edge
an involutive bijection.  For a state `x`, let `o(x)` and `e(x)` be its odd
and nonempty even admissible-cut counts.

**Target 22F3 (componentwise BK average).**  In the nontrivial sector
`dim O_p>P_p`, every connected component `C` of `K(d)` satisfies

```text
sum_(x in C)[o(x)-indicator_(o(x)>0)]
 <=sum_(x in C)e(x).                                  (P22.30)
```

Summing `(P22.30)` over the components gives Target 22F2.  Unlike the false
fixed-`alpha` statement, Target 22F3 allows exactly the cross-`alpha`
transfers generated by adjacent tableau switching; unlike an unrestricted
average, it presents a local graph on which a spanning forest or discrete
Morse matching could prove the inequality.  The strict component enumerator
passed every source-deficit instance occurring among the `1,001`
support-disjoint cases with at most five factors, labels and targets at most
four, and total degree at most fourteen.

There is a useful exact matroid formulation of what a pairwise proof would
have to select.  For every state `x` with odd-cut set `O(x)`, form the
complete graph on `O(x)`.  A spanning tree has `o(x)-1` edges, exactly the
relation contribution of that fibre in `(P22.23)`.  Take the direct sum of
these graphic matroids over all states and orders.  On the other side, make
one resource for every nonempty even object `(x,S)`.  Any specified rule
which declares which resources may pay an edge gives a transversal matroid
on the same edge set.

**Lemma 22F4 (graphic--transversal certificate).**  If the graphic and
transversal matroids have a common independent set of size

```text
R_rel=sum_(prec,x)(o(x)-indicator_(o(x)>0)),          (P22.30A)
```

then Target 22F2 holds.

**Proof.**  The rank of the graphic summand at `x` is `o(x)-1` when
`o(x)>0` and zero otherwise.  Equality with the sum of all these ranks
forces the selected edges to contain a spanning tree in every nontrivial
odd fibre.  Transversal independence matches the selected edges to distinct
positive resources.  Hence `R_rel` is at most the total number of positive
resources, which is the right side of `(P22.24)`.  QED.

This formulation sharply separates the missing issue from the already false
all-pairs bound.  Unfortunately, even a generous pair-local resource rule is
still too narrow.  Let an edge `{A,B}` use either any positive object in its
own fibre or any positive channel which is a parity-even union of the four
Boolean atoms

```text
A intersect B,  A\B,  B\A,  and  (A union B)^c.       (P22.30B)
```

The strict augmenting-path matroid-intersection enumerator finds the first
support-disjoint failure at

```text
remainder=[1^-,1^-,2^+,2^+,3^+],       target=1.      (P22.30C)
```

On the signed-class quotient there are `180` orders.  The odd source has
dimension eight per order, the rank sum is `876`, and therefore the exact
relation demand is `180*8-876=564`.  The positive supply is `180*5=900`, so
Target 22F2 has slack `336`.  Nevertheless XOR-only spanning trees have
maximum common rank `484`, and allowing both local objects and every
parity-even Boolean resolution in `(P22.30B)` raises it only to `548<564`.
Even replacing the two-cut Boolean algebra by the algebra of membership
atoms generated by the *entire odd fibre* leaves the same maximum `548`.
The full physical BK graph for this case is connected and has `2,256`
relations, `3,600` positive objects, and slack `1,344`, exactly four times
the quotient counts.  Thus the surviving component inequality cannot be
proved by a differential whose destination depends only on one colliding
pair of cuts and its Boolean regions; merely seeing all Boolean membership
atoms of the fibre is still insufficient.  It must use a non-Boolean
higher syzygy or genuine transport through other BK states.

The component picture itself is nontrivial.  For degree vector `[2,2,2,2]`,
`K(d)` has components of sizes `48` and `24`; only the first contains a
nontrivial admissible cut.  For `[2,2,3,3]` the same sizes occur, showing that
the gcd of the row multiplicities is not a complete component invariant.
For the first all-degree-at-least-two source deficit,

```text
remainder=[2^+,2^+,3^-,4^+,5^+],       target=2,
```

the two component sizes are `15,120` and `720`, and again only the first is
loaded.  The strict `gt-bk-single-loaded` sweep passed all `1,868`
support-disjoint source-deficit cases with at most six factors, labels and
targets at most four, and total degree at most sixteen.  This property is
special to the relevant signed sector: the unsigned degree vector
`[2,2,2,2,2,2]` already has two different components containing admissible
cuts.  Hence a proof may use the source-deficit and support-disjoint sign
hypotheses, but cannot assume a universal unique decomposable orbit.

The fibres also admit a matching-theoretic model which exposes where a
non-Boolean syzygy could live.  For a state `(prec,alpha)`, make
`alpha_v` top tokens and `beta_v=d_v-alpha_v` bottom tokens at every vertex.
Join a top token at `u` to a bottom token at `v` precisely when `u` precedes
`v` in `prec`.  This is a Ferrers bipartite graph `H_(prec,alpha)`.

**Lemma 22F5 (Ferrers separation form).**  A cut `T` is admissible for
`(prec,alpha)` if and only if both vertex-induced bipartite graphs

```text
H[T]  and  H[T^c]
```

have perfect matchings, where the target belongs to `T^c`.

**Proof.**  List the top and bottom tokens of one induced graph in
nondecreasing vertex order.  Because its adjacency matrix is Ferrers, it has
a perfect matching if and only if the `j`-th top token precedes the `j`-th
bottom token for every `j`; when this holds, pairing equal ranks is already
a perfect matching.  Equality of the two token counts is the terminal queue
condition, and the rank-by-rank inequalities are exactly the prefix queue
conditions in `(P22.16)`.  Apply this separately to `T` and `T^c`.  QED.

Thus Target 22F3 is an averaged signed uncrossing theorem for balanced
vertex separations of Ferrers graphs which remain perfectly matchable on
both sides.  The failure of the whole-fibre Boolean rule above says that the
standard lattice of membership atoms alone does not prove the required
uncrossing; a successful differential must also use alternating matching
paths, or equivalently the adjacent state transport `(P22.18)`.

The queue formulation has an exact finite-level analogue.  For an ordered
list `d_1,...,d_m` and a level-`k` fusion path

```text
h_0=0,h_1,...,h_m=0,
```

put

```text
alpha_i=(d_i+h_i-h_(i-1))/2,       beta_i=d_i-alpha_i. (P22.31)
```

**Lemma 22G (two-sided level queue).**  Level-`k` invariant fusion paths are
in bijection with integer vectors in `(P22.31)` satisfying, at every step,

```text
0<=alpha_i,beta_i<=d_i,
beta_i<=h_(i-1),                 alpha_i<=k-h_(i-1),
h_i=h_(i-1)-beta_i+alpha_i.                       (P22.32)
```

Thus `beta_i` consumes the lower queue while `alpha_i` consumes the upper
hole queue.  Omitting the upper inequality recovers the ordinary queue in
Lemma 22F.

**Proof.**  The ordinary Clebsch--Gordan conditions for
`V_(h_(i-1)) tensor V_(d_i)` are equivalent to the integrality in `(P22.31)`,
`0<=alpha_i,beta_i<=d_i`, and `beta_i<=h_(i-1)`.  The extra affine-wall
condition is

```text
h_(i-1)+d_i+h_i<=2k.
```

Substituting `h_i=h_(i-1)+d_i-2beta_i` turns it into
`h_(i-1)+alpha_i<=k`, which is the upper inequality in `(P22.32)`.
Every transformation is reversible.  QED.

Lemma 22G identifies the uniform finite target: construct the required
odd-to-positive matching on pairs of two-sided queue paths using moves that
preserve both inequalities in `(P22.32)`.  Such a matching would prove
`(PCP_k)` simultaneously for every `k`; dropping the upper queue would give
the ordinary matching studied in Targets 22C--22F.

Thus the remaining proof cannot be an order-independent polymatroid
exchange, an order-independent pairwise XOR injection, or a fibrewise
parity theorem.  The surviving statement is the adaptive cumulative Hall
matching of Target 22C2, or the weaker rank-plateau statement Target 22F1.

Two still stronger fixed-degeneration shortcuts are also false.  In the squarefree
Stanley--Reisner degeneration of the Plucker ring, crossing products vanish.
For remainder labels `[1,1,1,1]`, target `2`, and minus mask `3`, the odd
source dimension is four, its degenerate rank is one, and the nonempty even
supply is two.  Thus the special-fibre kernel has dimension three and is too
large by one.  In the domain Gelfand--Tsetlin degeneration a fixed vertex
order can also lose too much rank: for support-disjoint labels
`[1,1,2,2,3]`, target `1`, and minus mask `3`, the odd source dimension,
degenerate rank, and nonempty even supply are respectively `8,2,5`.
Reordering the six vertices as `[0,2,4,1,3,5]` raises the rank to six and
makes every tableau fibre satisfy `(P22.9)`.  Hence neither standard fixed
degeneration proves even the aggregate estimate automatically.

The strict exact C++ diagnostic
`character_ring_iter/analyze_su2_plucker_sr_parity.cpp` constructs both
degenerations and the tableau fibres.  Its order-adapted support-disjoint
fibre-order mode passed all 505 cases with at most six factors, remainder labels at most
three, targets at most three, and total degree at most fourteen.  A wider
replay passed all 2,422 support-disjoint cases with at most seven factors,
labels and targets at most four, and total degree at most sixteen.  Proposition
22B first occurs beyond that replay's total-degree cutoff.  Outside the
support-disjoint sector the fibrewise claim itself can fail for every vertex
order: six label-two factors, target two, and three minus copies have a fibre
with six odd and four even cuts.  Proposition 5 removes precisely this type
of mixed-sign equal-label case before Target 22B is invoked.
The earlier bounded replay is recorded in
`certificates/su2_plucker_gt_sd_best_4_7_4_16.log`.

The divided-difference factorization used earlier for the all-minus search
has a useful categorical and finite-level form.  In the ordinary
representation ring, or in `R_k` with `q<=k`, put

```text
Delta=e_1 tensor e_0-e_0 tensor e_1,
S_q=e_q tensor e_0+e_0 tensor e_q,
D_q=sum_(j=0)^(q-1) e_(q-1-j) tensor e_j.               (P23.1)
```

**Proposition 23 (fundamental-difference normal form).**  One has

```text
e_q tensor e_0-e_0 tensor e_q=Delta D_q.                (P23.2)
```

Consequently every signed remainder with `m` minus factors has the exact
form

```text
R=Delta^m H,
H=product_(i minus)D_(q_i) product_(i plus)S_(q_i).      (P23.3)
```

Here `H` is exchange-symmetric and belongs to the nonnegative tensor-product
basis cone.  Thus `(PCP)` and `(PCP_k)` reduce to the single rank-one target

```text
(id tensor tau)(Delta^m H) belongs to the nonnegative cone              (FDC)
```

for `H` in the multiplicative cone generated by the `D_q` and `S_q`.

**Proof.**  Write `e_q=U_q(e_1/2)`.  The divided-difference identity

```text
[U_q(X/2)-U_q(Y/2)]/(X-Y)
   =sum_(j=0)^(q-1)U_(q-1-j)(X/2)U_j(Y/2)               (P23.4)
```

follows by telescoping the recurrence
`XU_j(X/2)=U_(j+1)(X/2)+U_(j-1)(X/2)`.  Substituting
`X=e_1 tensor e_0` and `Y=e_0 tensor e_1` proves `(P23.2)`.
All indices in `D_q` are at most `q-1`; hence the same telescoping uses no
affine wall when `q<=k`, and the identity remains exact in every `R_k`.
Multiplying `(P23.2)` over the minus factors gives `(P23.3)`.  Every `D_q`
and `S_q` is symmetric and coefficientwise nonnegative, proving the final
assertions.  QED.

After the principal specialization used in Proposition 18, `(P23.2)` reads

```text
P_q(t)-t^(q/2)chi_q(y)
 =[1+t-t^(1/2)chi_1(y)]
   sum_(j=0)^(q-1)t^(j/2)P_(q-1-j)(t)chi_j(y).           (P23.5)
```

There is an exact supermodule behind both factors in `(P23.5)`.  Put
`a=t^(1/2)`, let `V=V_1`, and grade

```text
E=exterior^bullet(aV)=1 direct-sum aV direct-sum a^2 1.  (P23.6)
```

Here the displayed expression is the ordinary graded character; its
supercharacter is

```text
B(a,y)=1-a chi_1(y)+a^2.                                (P23.7)
```

**Lemma 23A (exterior--symmetric-power form).**  For every `q>=1`,

```text
Q_q(a,y)=sum_(j=0)^(q-1) a^j P_(q-1-j)(a^2) chi_j(y)
        =character Sym^(q-1)(E),                         (P23.8)

P_q(a^2)-a^q chi_q(y)=B(a,y)Q_q(a,y).                   (P23.9)
```

Consequently the outer polynomial `A_R` is the invariant super-Hilbert
polynomial of the genuine graded supermodule

```text
K_R=
 tensor_(i minus)[exterior^bullet(aV) tensor Sym^(q_i-1)(E)]
 tensor_(i plus)[(direct-sum_(j=0)^q_i a^(2j)1)
                  direct-sum a^(q_i)V_(q_i)],           (P23.10)
```

where parity occurs only in the first exterior factor of each minus
position.

**Proof.**  If `z,z^(-1)` are the weights of `V`, then

```text
sum_(n>=0) character Sym^n(E) u^n
 =1/[(1-u)(1-a^2u)(1-azu)(1-az^(-1)u)].                 (P23.11)
```

The first two denominator factors have coefficient `P_l(a^2)` in degree
`l`, and the last two have coefficient `a^j chi_j(y)` in degree `j`.
Convolution gives `(P23.8)`, while `(P23.9)` is `(P23.5)`.  Multiplying the
supercharacters of the factors in `(P23.10)` and extracting `SU(2)`
invariants gives exactly `(P18.3)`.  QED.

This supplies a concrete Koszul/Dirac object, but superdimension alone does
not give the required sign.  In particular, a natural attempt to use the
positive product formula for a higher-multiplicity rank-two Jacobi system
fails exactly.

For `m>=1`, put

```text
dnu_m=(x-y)^(2m)dmu(x)dmu(y).                            (P23.12)
```

Split `2m` minus factors into two groups of `m`.  The signed integral is

```text
< product_(i=1)^m D_(q_i),
  product_(plus j)S_(r_j) product_(i=m+1)^(2m)D_(q_i) >_(nu_m). (P23.13)
```

The measure `(P23.12)` is a rank-two `C_2/BC_2` Jacobi weight.  One could
therefore hope that every `m`-fold row product belongs to the nonnegative
cone of its triangular orthogonal basis and that every `S_q` preserves that
cone.  Both necessary statements are false at the first available degrees.

**Proposition 23B (`BC_2` hypergroup-cone obstruction).**  Let
`P_lambda^(m)` be obtained by Gram--Schmidt orthogonalization of the
`USp(4)` characters `Psi_lambda` against `(P23.12)`, with leading
coefficient one and the lower dominance terms taken first.

1.  At `m=1`, this is the ordinary `USp(4)` character basis, and

```text
S_2=Psi_(2,0)-Psi_(1,1)+Psi_(0,0).                      (P23.14)
```

Thus `S_2` has coefficient `-1` and does not preserve its nonnegative cone.

2.  At `m=2`,

```text
P_(1,1)^(2)=Psi_(1,1)+3/5,
<P_(1,1)^(2),P_(1,1)^(2)>_(nu_2)=42/5,                 (P23.15)

<D_1 D_3,P_(1,1)^(2)>_(nu_2)=-14/5.                   (P23.16)
```

Hence the coefficient of `P_(1,1)^(2)` in the two-row half-product
`D_1D_3=Psi_(2,0)` is `-1/3`.

**Proof.**  Formula `(P23.14)` is `(P8.2)` at `q=2`.  For the second
statement use `Psi_(1,1)=xy+1`, `Psi_(2,0)=x^2+xy+y^2-2`, and the semicircle
moments

```text
integral x^2 dmu=1,  integral x^4 dmu=2,
integral x^6 dmu=5.                                  (P23.17)
```

Direct expansion of `(x-y)^4` gives

```text
integral (x-y)^4 dmu^2=10,
integral (x-y)^4 Psi_(1,1)dmu^2=-6,
```

which yields the first line of `(P23.15)`.  The same expansion gives the
norm and pairing in `(P23.15)--(P23.16)`; their quotient is `-1/3`.
QED.

The strict exact rational C++ diagnostic
`character_ring_iter/analyze_su2_bc2_hypergroup.cpp` constructs the
triangular basis directly from the semicircle moments.  The commands with
arguments `1 4` and `2 4` reproduce respectively `(P23.14)` and
`(P23.16)`.  Proposition 23B does not contradict `Q3`: it rules out the
stronger proof in which generic Jacobi-hypergroup positivity is applied
after the half-product split.  Any successful use of `(P23.10)` or
`(P23.13)` must retain cancellations between several orthogonal types.

Thus all label dependence has moved into an explicitly positive quotient;
all cancellation is carried by the single fundamental factor `Delta^m`.
This is the natural place to seek the missing Temperley--Lieb differential
in Proposition 22.

In the ordinary category the positive quotient itself is the character of
one familiar module.  If

```text
A=V_1 tensor 1,                 B=1 tensor V_1,
```

then the symmetric-power decomposition gives

```text
D_q=character Sym^(q-1)(A direct-sum B).                 (P23.18)
```

Thus `(FDC)` is not asking for positivity of an artificial quotient: after
extracting `Delta^m`, every remaining minus factor is an actual polynomial
`SU(2) x SU(2)` module.  The unresolved sign is precisely the Euler index
created by the common virtual factor `(A-B)^m`.

The restriction on `H` in `(FDC)` is essential.  Positivity fails on the
whole symmetric nonnegative tensor cone: with

```text
m=1,                 H=e_1 tensor e_1,
```

one gets `(id tensor tau)(Delta H)=-e_1`.  The strict C++ diagnostic
`character_ring_iter/search_su2_fundamental_difference_cone.cpp` finds this
as its first symmetric-atom counterexample.  A proof must therefore preserve
the complete `D_q,S_q` generator groupings, just as the outer-Turan atom
counterexamples require complete irreducible strings.

For one minus factor, the normal form has an exact source--sink transfer
description.  Expand an arbitrary tensor character as

```text
H=sum_(b>=0) H_b tensor e_b,       H_(-1)=0,
L_b(H)=e_1 H_b-H_(b-1)-H_(b+1).                       (P24.1)
```

**Proposition 24 (fundamental transfer recurrence).**  One has

```text
(id tensor tau)(Delta H)=L_0(H).                       (P24.2)
```

Moreover, if `S_q=e_q tensor e_0+e_0 tensor e_q`, then

```text
L_b(H S_q)=e_q L_b(H)+sum_c N_(q,c)^b L_c(H),          (P24.3)
```

where `N_(q,c)^b` is the `SU(2)` tensor multiplicity.  The divided-difference
generator has the two-point boundary vector

```text
L_b(D_q)= e_q indicator_(b=0)-e_0 indicator_(b=q).      (P24.4)
```

Thus the one-minus case of `(PCP)` is exactly the assertion that positive
transfer by the operators

```text
(T_q L)_b=e_q L_b+sum_c N_(q,c)^b L_c                 (P24.5)
```

always makes the source at `b=0` dominate the sink at `b=q` in the zeroth
row.  The identical recurrence holds on the finite `A_(k+1)` fusion graph.

**Proof.**  Since `tau(e_1e_b)=indicator_(b=1)`, contraction of
`Delta H` gives `e_1H_0-H_1=L_0(H)`.  If `K=HS_q`, then

```text
K_b=e_qH_b+sum_c N_(q,c)^b H_c.
```

Apply `(P24.1)` and use associativity of the fusion coefficients, or
equivalently commutation of the fundamental adjacency matrix with every
fusion matrix `N_q`, to obtain `(P24.3)`.  Finally the `b`th row of `D_q` is
`e_(q-1-b)` for `0<=b<q` and zero otherwise.  The fundamental fusion
recurrence cancels every interior term, leaving exactly `(P24.4)`.  The
finite path graph obeys the same calculation because all indices used in
`D_q` lie below the affine wall.  QED.

Proposition 24 recovers the previously isolated two-minus boundary
domination problem, but now as a positive transfer of one explicit
source--sink pair.  A path injection proving that the source contribution
dominates the sink would settle the complete one-minus partial-character
sector simultaneously in ordinary `SU(2)` and all `SU(2)_k`.

The first previously open edge beyond Proposition 12E1 can now be closed
uniformly outside three bottom packet labels.  The point is that in the
support-disjoint problem the complete negative allocation load has only
three types, while either of the two types involving a channel `1` forces a
quadratic reservoir.

For a submultiset `B` of `L={q,r,s,t}`, write `m_B(c)` for the multiplicity
of `V_c` in its tensor product and put

```text
M=m_(p-1,q,r,s,t)(a).
```

Assume throughout this paragraph that

```text
2<=p<=q<=r<=s<=t,          a notin {p,q,r,s,t}.          (P25.1)
```

In the exact allocation formula `(P12E1.2)`, now with four suffix labels,
all singleton negative terms vanish.  Disjoint support kills the first
three-label negative term and both four-label negative terms.  Hence the
whole negative load is

```text
N=N_A+N_B+N_C,                                      (P25.2)

N_A=sum_(B subset L, |B|=2) m_B(p-1)m_(L minus B)(a),

N_B=sum_(B subset L, |B|=2) m_B(1)
                              m_((L minus B) union {p})(a),

N_C=sum_(B subset L, |B|=3) m_B(1)
                              indicator_(a in (L minus B) star p).
```

The empty allocation contributes `M`, and every omitted term is positive.
Thus `M>=N` is sufficient for the ordered packet coefficient.

**Lemma 25 (two tent-overlap bounds).**  Let `p>=2`.

1.  If `u,v,x,y>=p`, `|u-v|=1`, `a` is different from `x,y`, and
    `m_(x,y,p)(a)>0`, then

```text
m_(p-1,u,v,x,y)(a)>=p(p+2).                         (P25.3)
```

2.  If `u,v,w,x>=p`, `m_(u,v,w)(1)>0`, `a!=x`, and
    `a in x star p`, then

```text
m_(p-1,u,v,w,x)(a)>=p(p+2).                         (P25.4)
```

**Proof.**  We use only the interval form of the Clebsch--Gordan rule.  For
three labels set

```text
T_(b,c,d)(n)=#(I(b,c) intersection I(d,n)).          (P25.5)
```

For `(P25.3)`, the adjacent pair satisfies

```text
chi_u chi_v >= chi_p chi_(p+1)                       (P25.6)
```

coefficientwise.  It is therefore enough to bound

```text
<chi_(p-1)chi_p chi_(p+1), chi_a chi_x chi_y>.       (P25.7)
```

The interval `I(x,y)` has at least `p+1` points and, by the hypothesis in
`(P25.3)`, meets `I(a,p)`.  Choose `p+1` consecutive points of `I(x,y)`
containing a point of this intersection.  They form `I(c,p)` for some
`c>=p`, so `(P25.7)` is bounded below by

```text
<chi_(p-1)chi_p chi_(p+1), chi_a chi_c chi_p>.       (P25.8)
```

Put `d=(c-a)/2` and `h=(c+a)/2`.  The intersection condition says
`|d|<=p`; positivity of `a,c` gives `h-d>=1` and `h+d>=p`.  Pair the six
factors in `(P25.8)` as `(p-1,p+1)`, `(p,p)`, and `(a,c)`, and divide all
three intermediate labels by two.  The resulting invariant paths are

```text
1<=i<=p,       0<=j<=p,       |d|<=k<=h,
|i-j|<=k<=i+j.                                      (P25.9)
```

If `f_p(k)` denotes the number of `(i,j)` in the first two ranges satisfying
the last inequalities, then for `0<=k<=p+1`

```text
f_p(k)=p(p+1)-max(p-k,0)^2-k(k-1)/2.                 (P25.10)
```

This follows by deleting the disjoint regions `|i-j|>k` and `i+j<k` from
the `p` by `p+1` rectangle.  Summing `(P25.10)` under the three constraints
on `d,h` shows that the minimum is attained at `d=p,h=p+1`.  Terms with
`k>p+1` are nonnegative and cannot lower the sum.  At the minimum,

```text
f_p(p)+f_p(p+1)
 =p(p+3)/2+p(p+1)/2=p(p+2).                         (P25.11)
```

This proves `(P25.3)`.

For `(P25.4)`, choose `p+1` consecutive channels in `I(u,v)` containing
one of the channels in `I(w,1)`.  They again form `I(c,p)` with `c>=p`.
Moreover `p in I(a,x)`.  Since `a!=x`, `x>=p`, and the endpoints have the
same parity, the upper endpoint of `I(a,x)` is at least `p+2`; hence

```text
chi_a chi_x >= chi_p+chi_(p+2)=chi_1 chi_(p+1).       (P25.12)
```

It is therefore enough to count

```text
<chi_c chi_p chi_w, chi_(p-1)chi_1 chi_(p+1)>.       (P25.13)
```

Pair the factors as `(c,w)`, `(p-1,p+1)`, and `(p,1)`.  The last two
fusion intervals are

```text
{2,4,...,2p},                 {p-1,p+1}.              (P25.14)
```

The common-channel hypothesis says that the lower endpoint of `I(c,w)` is
at most `p+1`, with the required parity.  Among all such intervals with
`c,w>=p`, the number of invariant triples with the two intervals in
`(P25.14)` is minimized when the lower endpoint is `p+1`.  In that extremal
case `I(c,w)` contains every label from `p+1` through `3p+1`.  For the
channel `2i`, `1<=i<=p`, the choices paired with `p-1` and `p+1` number
respectively `i` and `i+1`.  Therefore the minimum is

```text
sum_(i=1)^p (i+i+1)=p(p+2).                          (P25.15)
```

Moving the lower endpoint inward only adds channels inside one of these
triangle intervals, so it cannot decrease the count.  This proves
`(P25.4)`.  QED.

There is also a sharp linear bound for `(P25.2)`.

**Lemma 25A (negative-load bound).**  Under `(P25.1)`,

```text
N<=4p+12.                                             (P25.16)
```

**Proof.**  Every term in `(P25.2)` is a quantifier-free Presburger
function.  In particular,

```text
m_(x,y,p)(a)=#(I(x,y) intersection I(a,p))            (P25.17)
```

is the cardinality of the intersection of two parity intervals.  The
strict C++/Z3 verifier
`character_ring_iter/verify_su2_opd_four_suffix_z3.cpp`, in mode
`--negative-bound`, encodes `(P25.1)`, all six two-label allocations, all
four three-label allocations, and the exact interval cardinality
`(P25.17)`.  It asks for an integer solution of `N>4p+12`.  Z3 4.8.12
returns `unsat`.  No label cutoff or bounded enumeration occurs in this
query.  QED.

**Proposition 25B (support-disjoint four-suffix edge for `p>=5`).**  If
`(P25.1)` holds and `p>=5`, then

```text
coefficient_(e_a wedge e_0)
 dGamma(N_t)dGamma(N_s)dGamma(N_r)dGamma(N_q)Y_p >=0. (P25.18)
```

**Proof.**  If `N_B+N_C>0`, one of the hypotheses of Lemma 25 occurs, so

```text
M>=p(p+2)>=4p+12>=N,                                 (P25.19)
```

where the middle inequality is equivalent to `(p-5)(p+3)>=0`.  If
`N_B=N_C=0`, every term of `N=N_A` is a distinct admissible `3|3` split of
the six labels `(p-1,q,r,s,t,a)`.  Proposition 1C gives `M>=N_A`.  In both
cases the empty positive allocation pays the complete negative load; all
other allocations are nonnegative surplus.  QED.

It remains to close the three fixed values `p=2,3,4`.  This can be done by
a second cutoff-free Presburger certificate.  The certificate is much
smaller after Lemma 25: it only needs to consider

```text
N>p(p+2).                                             (P25.20)
```

Indeed, if `(P25.20)` fails and `N_B+N_C>0`, Lemma 25 gives `M>=N`; if
`N_B=N_C=0`, Proposition 1C does so.

Here is the explicit path pool used in the remaining chamber.  The first
channel in the fusion tree for `M` is

```text
x_i=q-p+1+2i,                         0<=i<p.          (P25.21)
```

Let `delta(s,t,a)` be the smallest label in the support of
`chi_s chi_t chi_a`, namely

```text
delta=max(2 max(s,t,a)-s-t-a, (s+t+a) mod 2).         (P25.22)
```

For each `x_i`, start at

```text
l_i=max(|x_i-r|,delta),
y_(i,j)=l_i+2j.                                       (P25.23)
```

Retain the first `h_p` valid values, where

```text
h_2=10,                 h_3=8,                 h_4=7. (P25.24)
```

For each retained `y_(i,j)`, put

```text
z_(i,j,k)=max(|y_(i,j)-s|,|a-t|)+2k,       0<=k<4,   (P25.25)
```

and retain it exactly when it lies below both upper fusion endpoints.
Parity agreement and the two upper bounds are linear predicates.  Every
retained triple `(x_i,y_(i,j),z_(i,j,k))` is a distinct fusion path counted
by `M`.

For `p=3,4`, let `C_p` be the number of retained paths.  The verifier asks
for a solution of

```text
(P25.1), (P25.20),                    C_p<N.          (P25.26)
```

For `p=2`, the selected positive reservoir also includes

```text
S=sum_(u in L, u=2) m_((L minus {u}) union {1})(a),
E=indicator_(a=1)[m_L(0)+m_L(2)].                    (P25.27)
```

The corresponding query is `C_2+S+E<N`.  Interval cardinalities in `S,E`
are encoded exactly; truncation merely selects genuine paths and is never
used as an equality for the full multiplicity.

The verifier partitions `p=2` into the endpoint, strict-minimum interior,
and repeated-minimum interior chambers.  It partitions `p=3,4` according
to `q=p` or `q>p`; the source now refines these further by endpoint and the
position of `a`, only to reduce solver time.  Every chamber returned
`unsat`:

```text
p=2: endpoint, strict interior, repeated-minimum interior;
p=3: q=p, q>p;
p=4: q=p, q>p.                                       (P25.28)
```

These are unbounded integer queries.  The constants in `(P25.24)` and
`(P25.25)` bound only the displayed witness pool; there is no bound on
`q,r,s,t,a`.

**Proposition 25C (complete support-disjoint four-suffix edge).**  Under
`(P25.1)`, the ordered packet coefficient `(P25.18)` is nonnegative for
every `p>=2`.

**Proof.**  Proposition 25B handles `p>=5`.  For `p=2,3,4`, first use
Lemma 25 or Proposition 1C unless `(P25.20)` holds.  In the remaining
chamber, the unsatisfiability statements `(P25.26)`--`(P25.28)` give
`C_p>=N` for `p=3,4` and `C_2+S+E>=N` for `p=2`.  These are counts of
distinct fusion paths inside the positive allocations, so the complete
packet coefficient is nonnegative.  QED.

This closes the first support-disjoint ordered packet case beyond
Proposition 12E1.  It is the four-suffix, seven-factor edge leaf.  General
ordered packet domination still requires a mechanism that survives
arbitrarily many suffix factors; Proposition 12E2 shows that neither a
cardinality prefix nor an outside-in shell filtration can provide it.

The strict exact C++ verifier
`character_ring_iter/verify_su2_torus_factorization.cpp` independently
compares `(P13.1)` with the subset invariant-multiplicity formula.  It passed
all 1,519 separately sorted signed multisets with labels at most four and at
most six factors.  This is a cross-check of the identity, not a finite proof
of the remaining inequality.

The analyzer mode also tested the auxiliary central-shell inequality

```text
B_(0,0)>=B_(2,2).
```

It passed all 63,090 separately sorted signed multisets with labels at most
six and at most eight factors.  This auxiliary statement remains unproved.
It cannot be strengthened to radial monotonicity: minus labels `[1,1]` and
plus labels `[1,1]` have `B_(2,2)=0<B_(4,4)=1`.  Nor does the shifted
three-charge inequality hold away from the bottom: plus labels `[1,1,4]`
give a negative shifted value at charge `2`.  Thus both Proposition 14 and
`(OTM)` are genuinely bottom/outer-boundary statements, rather than
translation-invariant inequalities on the whole charge lattice.

The individual feature terms in `(P13.1)` are not nonnegative.  With two
minus labels `1,3`, the feature choice `(j_1,j_2)=(1,3)` has

```text
A=(u-u^(-1))(u^3-u^(-3)),   q=-2,
```

while the choice `(1,1)` has `q=2`; their sum is zero, as required by Schur
orthogonality.  Hence the next exact target is a cumulative-weight theorem:
the nested sets `W_p` must compensate the negative atomic features.  Unlike
`(SAT)`, this formulation already includes arbitrary even minus parity and
reduces the Weyl obstruction to the charge levels `0,2,4`.
