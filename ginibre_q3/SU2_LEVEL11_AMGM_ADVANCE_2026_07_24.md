# Exact AM-GM program for `O_11`: completed

Date: 2026-07-24

This note is retained as the historical entry point for the weighted AM-GM
method. The program described here has now been completed: together with Hall,
braid, Farey-fan, translated-tail, cone-tree, separated-pole, and exact fusion
certificates, it exhausts all 27,962 residual regimes of `O_11`.

The full statement and exact census are in

```text
SU2_LEVEL11_FULL_GKS2_2026_07_24.md.
```

The principal weighted AM-GM lemma is unchanged. For

```text
F(r)=sum_(j in P) C_j product_l lambda_(j,l)^r_l
     -sum_(n in N) D_n product_l mu_(n,l)^r_l,
```

choose rational weights `alpha_(n,j)>=0` satisfying

```text
sum_j alpha_(n,j)=1,
product_j lambda_(j,l)^alpha_(n,j) >= mu_(n,l),
sum_n alpha_(n,j)D_n <= C_j.
```

Weighted AM-GM then gives `F(r)>=0` throughout the residual orthant. Exact
certificates use rational algebraic intervals for coefficients and capacities,
with 512-bit directed-rounding MPFR intervals for large logarithmic products.
Numerical optimization proposes rational tables but is not proof evidence.

The originally isolated chamber

```text
(B_1^-)^(2+2p)(B_5^+)^(1+2q)
```

remains independently replayable through

```text
character_ring_iter/verify_su2_o11_first_chamber_exact.cpp.
```

The complete source-and-ledger package was replayed before a container reset
but still needs to be restored to the repository archive. This reproducibility
status is recorded explicitly in the full theorem note.
