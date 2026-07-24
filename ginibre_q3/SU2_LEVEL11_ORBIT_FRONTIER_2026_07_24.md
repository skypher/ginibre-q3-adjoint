# `SU(2)_11` and the rank-six odd-orbit frontier

Date: 2026-07-24

## Exact reduction

By odd-level simple-current lifting,

```text
full SU(2)_11 scalar GKS2*  <=>  scalar GKS2* in the rank-six orbit ring O_11.
```

The full theorem remains open, but the former first obstruction is now
resolved exactly.

## First chamber proved

For all `p,q>=0`,

```text
[V_0 tensor V_0]
 (B_1 tensor 1-1 tensor B_1)^(2+2p)
 (B_5 tensor 1+1 tensor B_5)^(1+2q) >= 0.
```

Five denominator-100 weighted AM-GM certificates cover the infinite regions

```text
q>=4; q=3; q=2,p>=1; q=1,p>=1; q=0,p>=2,
```

and exact integral fusion arithmetic gives zero at the four omitted points.

## Further exact progress

The same rational-algebraic method proves 25 additional residual orthants by
26 exact certificates. One orthant is split into a finite strip and translated
tail. The strict four-thread replay reports

```text
SU2_O11_AMGM_FRONTIER_EXACT PASS
certificates=26 residual_keys=25 denominator=100 threads=4.
```

All roots are isolated by rational Sturm sequences, and every base and capacity
inequality is decided by rational interval endpoints. Numerical optimization
was used only to discover integer weight tables.

## Remaining target

Full `O_11`, and hence full `SU(2)_11`, remains unproved. The next work is to
classify the residual orthants not covered by Hall/fan transport or the current
AM-GM ledger, seeking translated-tail and endpoint-lattice certificates with
exact face checks.
