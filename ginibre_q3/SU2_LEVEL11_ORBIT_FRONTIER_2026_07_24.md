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
certificates and four exact zero leaves. The strict C++ proof is

```text
character_ring_iter/verify_su2_o11_first_chamber_exact.cpp.
```

## Further exact progress

The same weighted AM-GM mechanism proves:

```text
25 hand-audited residual keys, using 26 certificates;
1558 batch-ledger residual keys in 294 support/parity chambers.
```

Every accepted inequality is replayed with rational Sturm isolation and
rational interval arithmetic. The numerical allocation search is not proof
evidence.

The compact selector payload is

```text
certificates/su2_o11_amgm_batch_ledger.b64
```

with SHA-256

```text
4bad077b05a9f5554729da30af081728abe1578b70dc271a50813dab01794c69.
```

The proof structure and reproduction commands are in
`SU2_LEVEL11_AMGM_ADVANCE_2026_07_24.md`.

## Remaining frontier

A conservative reduced-budget transport census left 142 candidate regions
not covered by one global denominator-100 AM-GM allocation. This is not a
counterexample signal. The remaining regions may require lattice splitting,
translated tails with exact face certificates, a different transport, or a
structural Turan identity.

An authoritative full 100-ray/100-shift census is also still needed. Full
`O_11`, and hence full `SU(2)_11`, remains open.
