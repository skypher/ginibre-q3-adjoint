# Primary-source audit: adjacent nonabelian Ginibre inequalities

This record documents the source reading behind the related-work paragraph in
`full_q3_extension.tex`.  It deliberately records access and digests without
redistributing third-party PDFs in the publication source package.

## Garrett S. Sylvester (1980)

Primary source: G. S. Sylvester, *The Ginibre inequality*, Communications in
Mathematical Physics 73 (1980), no. 2, 105--114,
DOI `10.1007/BF01198120`.

Directly checked in the author-uploaded full text; its bibliographic identity
was cross-checked against the publisher DOI record:

1. Section 1, equation (1.1), formulates the Ginibre integral for a graph using
   two independent copies of sphere-valued spins and a sign on every edge.
2. Section 3 gives graph counterexamples in spin dimension at least three;
   Example 3.1 is the explicit analytic counterexample.
3. Section 4 and the appendix distinguish failure of the general Ginibre
   inequality from the still-unsettled Griffiths inequalities and report the
   numerical graph survey.  The appendix includes the exact negative value
   `-256/1366875` for graph 15.

The publisher metadata and author-uploaded full text were accessed on
2026-07-19 through `https://doi.org/10.1007/BF01198120` and
`https://www.researchgate.net/publication/38328963_The_Ginibre_inequality`.
The latter identifies the file as author content uploaded by Garrett
Sylvester.  The manuscript cites only the section-level conclusions above.

## Abdelmalek Abdesselam (2022)

Primary source: A. Abdesselam, *Non-Abelian correlation inequalities and
stable determinantal polynomials*, arXiv:2207.07603 (2022),
`https://arxiv.org/abs/2207.07603`.

The downloaded arXiv PDF has SHA-256
`ee535a5a320e11ff6553bd491ea54f40705bca15b44db66f4c8c3d165b197dea`.
The following loci were checked directly:

1. Section 1 states that invariant-observable GKS2 is known for `O(N)` at
   `N=1,2`, while the finite `N>=3` case remains open.  The same section
   defines the general Ginibre and padded general Ginibre families, recalls
   Sylvester's counterexamples, and states that they do not furnish a
   counterexample to the two-minus/CGKS2 subfamily.
2. Theorem 2.1 proves the large-power asymptotic for `O(N)` correlations, and
   Theorem 2.2 gives the corresponding projective-space asymptotic.
3. Theorem 2.3 proves all padded Ginibre inequalities for inverse
   half-integer powers of stable determinantal polynomials, used there as
   substitutes for the original finite-power correlations.

## Scope conclusion used in the paper

Sylvester has priority for nonabelian/higher-spin graph-cycle counterexamples.
The present paper's negative result is instead a particularly direct exact
witness for Ginibre's compact-group all-positive-definite candidate on
`SO(3)`.  Abdesselam's results are adjacent but concern large-power limits and
stable-determinantal substitutes.  Neither source proves the present
finite-rank adjoint-character Haar hierarchy, and the present theorem does not
claim to settle finite-`N` invariant-observable GKS2.

## Ilya Chevyrev and Christophe Garban (2025)

Primary source: I. Chevyrev and C. Garban, *Villain Action in Lattice Gauge
Theory*, Journal of Statistical Physics 192 (2025), Paper 38,
DOI `10.1007/s10955-025-03420-1`.

The open-access publisher version was checked directly at
`https://doi.org/10.1007/s10955-025-03420-1`:

1. Theorem 1.5 proves total-variation convergence of carpet-graph actions to
   the Villain lattice-gauge action for compact connected structure groups;
   the discussion explicitly includes nonabelian examples such as `SU(3)`.
2. Corollary 1.6 proves Wilson-loop monotonicity for the abelian `U(1)`
   Villain lattice-gauge theory.
3. Remark 1.7 explains the relation of that abelian monotonicity to the
   Ginibre inequality, while the nonabelian part of the paper is the
   Villain-action limiting construction rather than a nonabelian Ginibre
   correlation theorem.

This source therefore supports the paper's application boundary: rigorous
nonabelian lattice-gauge Villain limits exist, but adjoint representation
content alone does not turn the present Haar theorem into a nonabelian
lattice-gauge correlation inequality.
