# `SU(2)_11` and the rank-six odd-orbit frontier

Date: 2026-07-24

## Exact reduction

By odd-level simple-current lifting,

```text
full SU(2)_11 scalar GKS2*  <=>  scalar GKS2* in the rank-six orbit ring O_11.
```

Once the scalar theorem is proved, adjoining one plus factor gives the
complete partial-character column. The remaining odd-level problem is
therefore entirely the rank-six orbit ring.

Use the even-lift basis

```text
B_0=V_0, B_1=V_2, B_2=V_4, B_3=V_6, B_4=V_8, B_5=V_10.
```

## First chamber resolved

The former first obstruction

```text
(B_1^-)^(2+2p) (B_5^+)^(1+2q),       p,q>=0,
```

is now proved for all exponents by five exact denominator-100 AM-GM
certificates and four exact zero leaves. The authoritative strict C++ proof is

```text
character_ring_iter/verify_su2_o11_first_chamber_exact.cpp.
```

## Supplementary exact progress

The same turn produced locally exact rational C++ replays for

```text
25 hand-audited residual keys, using 26 certificates;
1558 batch-generated residual keys in 294 support/parity chambers.
```

The candidate allocations were generated with four parallel MILP workers but
accepted only after rational Sturm isolation and rational interval replay. The
numerical allocation search itself is not proof evidence.

The supplementary selector payload was not added to the authoritative archive
in this turn because the connector could not reliably upload the compact
90 KiB object. Its run summary and hashes are recorded in
`certificates/su2_o11_amgm_exact.log`; the authoritative repository theorem is
the complete first-chamber proof above.

## Remaining frontier

A conservative reduced-budget transport census left 142 candidate regions
not covered by one global denominator-100 AM-GM allocation. This is not a
counterexample signal. The remaining regions may require lattice splitting,
translated tails with exact face certificates, a different transport, or a
structural Turan identity.

An authoritative full 100-ray/100-shift census is also still needed. Full
`O_11`, and hence full `SU(2)_11`, remains open.
